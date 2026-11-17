use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use eimdl_native::fasta;
use eimdl_native::ffi::{self, MosaicFlags, SyncmerIndex};
use eimdl_native::hgs2::{self, CompressOptions, MODE_COLLECTION, MODE_REF, MODE_SINGLE};

fn make_temp_dir(tag: &str) -> PathBuf {
    let nanos = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let p = std::env::temp_dir().join(format!("eimdl_{tag}_{}_{}", std::process::id(), nanos));
    fs::create_dir_all(&p).unwrap();
    p
}

fn lcg_dna(seed: u64, len: usize) -> Vec<u8> {
    let mut x = seed;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        out.push(((x >> 32) & 3) as u8);
    }
    out
}

#[test]
fn deflate_and_zlib_roundtrip() {
    let sample = b">chr1 synthetic test header\nACGTACGTACGTACGTNNNNRYKM\r\n>chr2\nacgtacgt\n".repeat(40);
    let comp = ffi::deflate_raw(&sample).unwrap();
    let dec = ffi::inflate_raw(&comp, sample.len()).unwrap();
    assert_eq!(dec, sample);

    let zc = ffi::zlib_compress(&sample).unwrap();
    let zd = ffi::zlib_decompress(&zc).unwrap();
    assert_eq!(zd, sample);
}

#[test]
fn novelty_crm_turbo_and_law4_bypass() {
    let mut rep = Vec::with_capacity(80_000);
    let motif = lcg_dna(0x11223344, 512);
    for _ in 0..150 {
        rep.extend_from_slice(&motif);
    }
    let comp = ffi::novelty_encode(&rep, 2, false).unwrap();
    let dec = ffi::novelty_decode(&comp, rep.len() as u64, 2, false).unwrap();
    assert_eq!(dec, rep);
    assert!(comp.len() * 4 < rep.len(), "Expected < 0.5 bpb on tandem repeat, got {} / {}", comp.len(), rep.len());

    // Law IV cryptographic noise: raw 2-bit safety ceiling
    let noise = lcg_dna(0xDEADBEEFCAFEBABE, 200_000);
    let comp_n = ffi::novelty_encode(&noise, 2, false).unwrap();
    let dec_n = ffi::novelty_decode(&comp_n, noise.len() as u64, 2, false).unwrap();
    assert_eq!(dec_n, noise);
    let bpb = (comp_n.len() as f64 * 8.0) / (noise.len() as f64);
    assert!(bpb <= 2.01, "Expected Law IV ceiling <= 2.01 bpb, got {bpb}");
}

#[test]
fn mosaic_panel_and_self_roundtrip() {
    let donor0 = lcg_dna(0xABCD1234, 8000);
    let mut donor1 = donor0.clone();
    // Inject Kimura transitions + micro-indel
    for i in (100..7500).step_by(200) {
        donor1[i] ^= 2;
    }
    let mut tgt = donor1[..2000].to_vec();
    tgt.extend_from_slice(&lcg_dna(0x998877, 15)); // 15bp insertion splice
    tgt.extend_from_slice(&donor1[2010..]);
    // Also duplicate a reverse-complement block for Mode S self-repeat testing
    let tail_rc: Vec<u8> = donor0[1000..1500].iter().rev().map(|&b| b ^ 3).collect();
    tgt.extend_from_slice(&tail_rc);

    let idx = SyncmerIndex::build(&[donor0.as_slice(), donor1.as_slice()], 0, true, 2);
    let enc = ffi::encode_mosaic(&tgt, Some(&idx), 2, 1, &[1, 0], MosaicFlags::default(), 2).unwrap();
    let dec = ffi::decode_mosaic(&enc.mosaic_bytes, &enc.novelty_bases, &[donor0.as_slice(), donor1.as_slice()]).unwrap();
    assert_eq!(dec, tgt);
    assert!(enc.panel_match_bases > 6000);
}

#[test]
fn hgs2_end_to_end_modes_s_r_c() {
    let tmp = make_temp_dir("hgs2");
    let raw_fa = tmp.join("pangenome.fa");
    let arc_s = tmp.join("mode_s.hgs2");
    let arc_c = tmp.join("mode_c.hgs2");
    let dec_s = tmp.join("out_s.fa");
    let dec_c = tmp.join("out_c.fa");

    let base_seq: String = lcg_dna(0x7777, 18_000).iter().map(|&b| b"ACGT"[b as usize] as char).collect();
    let mut fasta_txt = Vec::new();
    for g in 0..4 {
        fasta_txt.extend_from_slice(format!(">strain_{g} synthetic pangenome\n").as_bytes());
        let mut s: Vec<u8> = base_seq.as_bytes().to_vec();
        if g > 0 {
            s[500 + g * 111] = b'a'; // lowercase soft-mask
            s[600 + g * 10] = b'N';  // N run
            s[700 + g * 7] = b'R';   // IUPAC ambiguity
        }
        for chunk in s.chunks(80) {
            fasta_txt.extend_from_slice(chunk);
            fasta_txt.extend_from_slice(if g % 2 == 0 { b"\n" } else { b"\r\n" });
        }
    }
    fs::write(&raw_fa, &fasta_txt).unwrap();

    // Mode S
    let opts_s = CompressOptions {
        mode: MODE_SINGLE as i32,
        single_genome: true,
        ..Default::default()
    };
    let rep_s = hgs2::compress_file(&raw_fa, &arc_s, &opts_s).unwrap();
    assert_eq!(rep_s.mode, MODE_SINGLE);
    hgs2::decompress_file(&arc_s, &dec_s, None, 2, false).unwrap();
    assert_eq!(fs::read(&dec_s).unwrap(), fasta_txt);

    // Mode C
    let opts_c = CompressOptions {
        mode: MODE_COLLECTION as i32,
        bounded: true,
        ..Default::default()
    };
    let rep_c = hgs2::compress_file(&raw_fa, &arc_c, &opts_c).unwrap();
    assert_eq!(rep_c.mode, MODE_COLLECTION);
    hgs2::decompress_file(&arc_c, &dec_c, None, 2, true).unwrap();
    assert_eq!(fs::read(&dec_c).unwrap(), fasta_txt);

    // Mode R + tamper pinning rejection
    let pdir = tmp.join("panel");
    fs::create_dir_all(&pdir).unwrap();
    let mut ref_bytes = b">donor0\n".to_vec();
    ref_bytes.extend_from_slice(base_seq.as_bytes());
    ref_bytes.push(b'\n');
    fs::write(pdir.join("donor0.fa"), &ref_bytes).unwrap();
    let tgt_fa = tmp.join("tgt.fa");
    let mut tgt_bytes = b">target0\r\n".to_vec();
    tgt_bytes.extend_from_slice(base_seq.as_bytes());
    tgt_bytes.extend_from_slice(b"ACGTacgtNNNNR\r\n");
    fs::write(&tgt_fa, &tgt_bytes).unwrap();

    let arc_r = tmp.join("mode_r.hgs2");
    let dec_r = tmp.join("dec_r.fa");
    let opts_r = CompressOptions {
        reference_panel: Some(pdir.clone()),
        ..Default::default()
    };
    let rep_r = hgs2::compress_file(&tgt_fa, &arc_r, &opts_r).unwrap();
    assert_eq!(rep_r.mode, MODE_REF);
    hgs2::decompress_file(&arc_r, &dec_r, Some(&pdir), 2, false).unwrap();
    assert_eq!(fs::read(&dec_r).unwrap(), tgt_bytes);

    // Tamper panel byte -> must reject SHA-256 manifest!
    let mut tampered = ref_bytes.clone();
    tampered[12] = b'T';
    fs::write(pdir.join("donor0.fa"), &tampered).unwrap();
    assert!(hgs2::decompress_file(&arc_r, &dec_r, Some(&pdir), 2, false).is_err());

    let info = hgs2::inspect_file(&arc_c).unwrap();
    assert_eq!(info.raw_bytes, fasta_txt.len() as u64);

    let _ = fs::remove_dir_all(&tmp);
    let _ = fasta::partition_file; // keep module linkage verified
}
