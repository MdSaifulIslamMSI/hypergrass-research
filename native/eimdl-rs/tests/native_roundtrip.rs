use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use eimdl_native::compute;
use eimdl_native::fasta;
use eimdl_native::ffi::{self, MosaicFlags, SyncmerIndex};
use eimdl_native::hgs2::{self, CompressOptions, MODE_COLLECTION, MODE_REF, MODE_SINGLE};
use eimdl_native::hgs3::{self, Hgs3CompressOptions, OmicInputKind};

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

    for sz in [0usize, 1, 12, 13, 35, 36, 37, 64] {
        let mut raw = vec![0u8; sz];
        for (i, b) in raw.iter_mut().enumerate() {
            *b = ((i * 73 + 17) & 0xFF) as u8;
        }
        let zc_s = ffi::zlib_compress(&raw).unwrap();
        let zd_s = ffi::zlib_decompress(&zc_s).unwrap_or_else(|e| panic!("zlib fail on len {sz}: {e}"));
        assert_eq!(zd_s, raw);
    }
    let mut cds_like = vec![0u8, 0, 0, 244, 3];
    cds_like.extend_from_slice(&ffi::sha256(b""));
    let zc_cds = ffi::zlib_compress(&cds_like).unwrap();
    let zd_cds = ffi::zlib_decompress(&zc_cds).unwrap_or_else(|e| panic!("zlib fail on cds_like (len={}): {e}", zc_cds.len()));
    assert_eq!(zd_cds, cds_like);
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

#[test]
fn ans4_dna_and_qual_conditioning_with_law4() {
    let mut rep = Vec::with_capacity(96_000);
    let motif = lcg_dna(0x31415926, 384);
    for _ in 0..250 {
        rep.extend_from_slice(&motif);
    }
    let comp = ffi::ans4_encode_dna(&rep).unwrap();
    let dec = ffi::ans4_decode_dna(&comp, rep.len()).unwrap();
    assert_eq!(dec, rep);
    let rep_bpb = (comp.len() as f64 * 8.0) / (rep.len() as f64);
    assert!(rep_bpb < 1.2, "Expected 4x rANS Order-4/6 < 1.2 bpb on repeat, got {rep_bpb}");

    // Multi-segment parallel rANS novelty (> 1.1M bases) + thread determinism
    let mut big = Vec::with_capacity(1_150_000);
    for _ in 0..3000 {
        big.extend_from_slice(&motif);
    }
    let big_1t = ffi::ans_novelty_encode(&big, 1).unwrap();
    let big_2t = ffi::ans_novelty_encode(&big, 2).unwrap();
    assert_eq!(big_1t, big_2t, "ANS novelty payload must be thread-count deterministic");
    let big_dec = ffi::ans_novelty_decode(&big_2t, big.len() as u64, 2).unwrap();
    assert_eq!(big_dec, big);

    // Law IV cryptographic noise ceiling
    let noise = lcg_dna(0xFEEDFACECAFE0001, 120_000);
    let noise_c = ffi::ans_novelty_encode(&noise, 2).unwrap();
    let noise_d = ffi::ans_novelty_decode(&noise_c, noise.len() as u64, 2).unwrap();
    assert_eq!(noise_d, noise);
    let noise_bpb = (noise_c.len() as f64 * 8.0) / (noise.len() as f64);
    assert!(noise_bpb <= 2.02, "Expected Law IV ceiling <= 2.02 bpb, got {noise_bpb}");

    // Pangenome mismatch + homopolymer conditioned Phred quality coder (including escape values > 96)
    let homo = ffi::compute_homo_bins(&rep[..40_000]);
    let mut mm = vec![0u8; 40_000];
    let mut quals = vec![b'I'; 40_000];
    for (i, q) in quals.iter_mut().enumerate() {
        let cycle = i % 150;
        if i % 37 == 0 {
            mm[i] = 1;
            *q = b'#'; // Low quality correlates with sequencing mismatch!
        } else if homo[i] >= 2 {
            *q = b'5';
        } else if cycle > 135 {
            *q = b'?';
        } else if i % 997 == 0 {
            *q = b'~'; // > 96 escape byte!
        }
    }
    let qc = ffi::ans4_encode_qual(&quals, &mm, &homo, 150).unwrap();
    let qd = ffi::ans4_decode_qual(&qc, &mm, &homo, quals.len(), 150).unwrap();
    assert_eq!(qd, quals);
    assert!(qc.len() < quals.len() / 4, "Expected high Phred compression, got {} / {}", qc.len(), quals.len());
}

#[test]
fn hgs3_fastq_lossless_and_panel_pinning() {
    let tmp = make_temp_dir("hgs3_fq");
    let pdir = tmp.join("panel");
    fs::create_dir_all(&pdir).unwrap();

    let ref_dna = lcg_dna(0x42424242, 6000);
    let ref_ascii: Vec<u8> = ref_dna.iter().map(|&b| b"ACGT"[b as usize]).collect();
    let mut ref_fa = b">donor_ref\n".to_vec();
    ref_fa.extend_from_slice(&ref_ascii);
    ref_fa.push(b'\n');
    fs::write(pdir.join("ref0.fa"), &ref_fa).unwrap();

    // Synthesize byte-diverse FASTQ with SNPs + homopolymers + IUPAC 'N' + lowercase + mixed +/CRLF + missing trailing EOF
    let mut raw_fq = Vec::new();
    for r in 0..80 {
        let id = format!("SRR999001.{r} inst:1:1101:{}:{} length=150", 1000 + r, 2000 + r * 3);
        let mut read = ref_ascii[r * 50..r * 50 + 150].to_vec();
        let mut qual = vec![b'I'; 150];
        if r % 5 == 0 {
            read[30] = b'T';
            if ref_ascii[r * 50 + 30] == b'T' {
                read[30] = b'A';
            }
            qual[30] = b'!';
        }
        if r % 9 == 0 {
            read[10] = b'N';
            read[11] = b'a';
        }
        let eol = if r % 2 == 0 { b"\n".as_slice() } else { b"\r\n".as_slice() };
        raw_fq.extend_from_slice(b"@");
        raw_fq.extend_from_slice(id.as_bytes());
        raw_fq.extend_from_slice(eol);
        raw_fq.extend_from_slice(&read);
        raw_fq.extend_from_slice(eol);
        raw_fq.extend_from_slice(b"+");
        if r == 3 {
            raw_fq.extend_from_slice(id.as_bytes());
        } else if r == 7 {
            raw_fq.extend_from_slice(b"custom_comment");
        }
        raw_fq.extend_from_slice(eol);
        raw_fq.extend_from_slice(&qual);
        if r + 1 < 80 {
            raw_fq.extend_from_slice(eol);
        }
    }

    let fq_in = tmp.join("reads.fastq");
    let arc3 = tmp.join("reads.hgs3");
    let fq_out = tmp.join("decoded.fastq");
    fs::write(&fq_in, &raw_fq).unwrap();

    let mut opts = Hgs3CompressOptions {
        format: OmicInputKind::Auto,
        emit_cds_digest: true,
        ..Default::default()
    };
    opts.base.reference_panel = Some(pdir.clone());
    let rep = hgs3::compress_omics(&fq_in, &arc3, &opts).unwrap();
    assert!(rep.is_fastq);
    assert_eq!(rep.mode, MODE_REF);
    assert!(rep.archive_bytes < raw_fq.len() as u64 / 3);

    hgs3::decompress_omics(&arc3, &fq_out, Some(&pdir), 2, false).unwrap();
    assert_eq!(fs::read(&fq_out).unwrap(), raw_fq, "HGS3 FASTQ decompression must be 100% byte-exact");

    // Compressed-domain analyze must work natively on HGS3-FASTQ!
    let ana = compute::analyze_archive(&arc3, Some(&pdir), 2, false).unwrap();
    assert_eq!(ana.format, "HGS3-FASTQ");
    assert!(ana.stats.panel_match_bases > 8000);

    // Panel tamper pinning rejection
    let mut bad = ref_fa.clone();
    bad[15] = if bad[15] == b'A' { b'C' } else { b'A' };
    fs::write(pdir.join("ref0.fa"), &bad).unwrap();
    assert!(hgs3::decompress_omics(&arc3, &fq_out, Some(&pdir), 2, false).is_err());

    // Also test HGS3-FASTA + CDS proteome ORF detection
    let mut fa_cds = b">orf_test\nATG".to_vec();
    fa_cds.extend_from_slice(&b"GCT".repeat(40)); // 40 Alanines -> > 30 aa ORF
    fa_cds.extend_from_slice(b"TAA\n");
    let fa_in = tmp.join("cds.fa");
    let fa_arc3 = tmp.join("cds.hgs3");
    let fa_out = tmp.join("cds_dec.fa");
    fs::write(&fa_in, &fa_cds).unwrap();
    let opts_fa = Hgs3CompressOptions {
        format: OmicInputKind::Fasta,
        emit_cds_digest: true,
        base: CompressOptions { single_genome: true, mode: MODE_SINGLE as i32, ..Default::default() },
    };
    let rep_fa = hgs3::compress_omics(&fa_in, &fa_arc3, &opts_fa).unwrap();
    assert!(!rep_fa.is_fastq);
    assert!(rep_fa.cds.unwrap().orf_count >= 1);
    hgs3::decompress_omics(&fa_arc3, &fa_out, None, 2, false).unwrap();
    assert_eq!(fs::read(&fa_out).unwrap(), fa_cds);

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn compressed_domain_analyze_variants_kmer_and_query() {
    let tmp = make_temp_dir("cd_genomics");
    let pdir = tmp.join("panel");
    fs::create_dir_all(&pdir).unwrap();

    let donor0 = lcg_dna(0xAAAA1111, 10_000);
    let donor1 = lcg_dna(0xBBBB2222, 10_000);
    let to_ascii = |v: &[u8]| -> Vec<u8> { v.iter().map(|&b| b"ACGT"[b as usize]).collect() };

    let mut d0_fa = b">donor0\n".to_vec();
    d0_fa.extend_from_slice(&to_ascii(&donor0));
    d0_fa.push(b'\n');
    let mut d1_fa = b">donor1\n".to_vec();
    d1_fa.extend_from_slice(&to_ascii(&donor1));
    d1_fa.push(b'\n');
    fs::write(pdir.join("d0.fa"), &d0_fa).unwrap();
    fs::write(pdir.join("d1.fa"), &d1_fa).unwrap();

    // Build target Mosaic genome:
    // Segment A: donor0[0..3000] with 4 SNPs (mod-3 codon wobble spacing + transition/transversion)
    // Segment B: 8bp insertion literal
    // Segment C: donor0[3005..5500] (implies 5bp deletion splice!)
    // Segment D: reverse-complement donor1[2000..5000] with 2 SNPs
    let mut seg_a = donor0[0..3000].to_vec();
    seg_a[500] ^= 2; // transition
    seg_a[503] ^= 2; // mod-3 transition
    seg_a[506] ^= 1; // mod-3 transversion
    seg_a[1200] ^= 3; // transversion
    let ins_lit = vec![0u8, 1, 2, 3, 2, 1, 0, 3];
    let seg_c = donor0[3005..5500].to_vec();
    let mut seg_d: Vec<u8> = donor1[2000..5000].iter().rev().map(|&b| b ^ 3).collect();
    seg_d[800] ^= 2;
    seg_d[803] ^= 1;

    let mut tgt_2bit = Vec::new();
    tgt_2bit.extend_from_slice(&seg_a);
    tgt_2bit.extend_from_slice(&ins_lit);
    tgt_2bit.extend_from_slice(&seg_c);
    tgt_2bit.extend_from_slice(&seg_d);
    let tgt_ascii = to_ascii(&tgt_2bit);

    let mut tgt_fa = b">recomb_target\n".to_vec();
    tgt_fa.extend_from_slice(&tgt_ascii);
    tgt_fa.push(b'\n');
    let tgt_path = tmp.join("tgt.fa");
    fs::write(&tgt_path, &tgt_fa).unwrap();

    let arc2 = tmp.join("tgt.hgs2");
    let arc3 = tmp.join("tgt.hgs3");
    let opts = CompressOptions { reference_panel: Some(pdir.clone()), ..Default::default() };
    hgs2::compress_file(&tgt_path, &arc2, &opts).unwrap();
    hgs3::compress_omics(
        &tgt_path,
        &arc3,
        &Hgs3CompressOptions { format: OmicInputKind::Fasta, base: opts.clone(), emit_cds_digest: false },
    )
    .unwrap();

    for arc in [&arc2, &arc3] {
        let ana = compute::analyze_archive(arc, Some(&pdir), 2, false).unwrap();
        assert_eq!(ana.total_bases, tgt_2bit.len() as u64);
        assert!(ana.stats.panel_match_bases > 7000);
        assert!(ana.stats.snps >= 4, "Expected compressed-domain SNPs, got {}", ana.stats.snps);
        assert!(ana.stats.recombination_switches >= 1);

        let vcf_out = tmp.join("calls.vcf");
        let vcf = compute::extract_variants_vcf(arc, Some(&pdir), Some(&vcf_out), 2, false).unwrap();
        assert!(vcf.vcf_text.contains("##fileformat=VCFv4.2"));
        assert!(vcf.snvs >= 4);
        assert!(vcf.vcf_text.contains("STRAND=-"), "Expected reverse-complement strand VCF call");

        // Brute-force k-mer verification on conserved, SNP, splice-junction, and palindrome k-mers!
        let queries = [
            (&tgt_ascii[100..115], false),  // inside conserved donor0
            (&tgt_ascii[498..512], true),   // spanning multiple SNPs
            (&tgt_ascii[2995..3015], true), // crossing insertion splice boundary!
            (&tgt_ascii[6295..6312], true), // inside reverse-complement donor1 + SNP!
            (b"ACGTACGT", true),            // palindrome across whole genome
        ];
        for (q_bytes, both) in queries {
            let q_str = std::str::from_utf8(q_bytes).unwrap();
            let q_rc: Vec<u8> = q_bytes
                .iter()
                .rev()
                .map(|&b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                })
                .collect();
            let mut expected = Vec::new();
            for pos in 0..=tgt_ascii.len() - q_bytes.len() {
                if &tgt_ascii[pos..pos + q_bytes.len()] == q_bytes {
                    expected.push((pos as u64, 0u8));
                }
                if both && &tgt_ascii[pos..pos + q_bytes.len()] == q_rc.as_slice() {
                    expected.push((pos as u64, 1u8));
                }
            }
            expected.sort();
            expected.dedup();

            let rep = compute::search_kmers(arc, Some(&pdir), q_str, both, 2, false).unwrap();
            let got: Vec<(u64, u8)> = rep.hits.iter().map(|h| (h.tgt_pos, h.strand)).collect();
            assert_eq!(got, expected, "K-mer mismatch for {q_str} (both={both}) on {}", arc.display());
        }

        let sliced = compute::query_slice(arc, Some(&pdir), 0, 2990, 3030, 2, false).unwrap();
        assert_eq!(sliced.as_bytes(), &tgt_ascii[2990..3030]);
    }

    let _ = fs::remove_dir_all(&tmp);
}
