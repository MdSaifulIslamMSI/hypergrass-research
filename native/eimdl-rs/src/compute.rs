use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::ffi::{self, HgsMosaicStats, HgsSpanEvent, HgsVariantRecord};
use crate::hgs2::{
    self, MAGIC, MAGIC3, MODE_COLLECTION, MODE_REF, ParsedHgs2, parse_container,
    verify_and_load_panel_for_decode,
};
use crate::hgs3::{self, CdsSummary, deserialize_cds_summary, parse_hgs3_fastq};

pub struct UnitInspectData {
    pub unit_index: usize,
    pub stats: HgsMosaicStats,
    pub variants: Vec<HgsVariantRecord>,
    pub spans: Vec<HgsSpanEvent>,
    pub var_seq_pool: Vec<u8>,
}

pub struct AnalyzeReport {
    pub format: &'static str,
    pub is_fastq: bool,
    pub mode: u8,
    pub raw_bytes: u64,
    pub total_bases: u64,
    pub stats: HgsMosaicStats,
    pub donor_names: Vec<String>,
    pub donor_match_bases: Vec<u64>,
    pub unit_count: usize,
    pub cds: Option<CdsSummary>,
    pub seconds: f64,
}

impl AnalyzeReport {
    pub fn ti_tv(&self) -> f64 {
        if self.stats.transversions == 0 {
            if self.stats.transitions > 0 { f64::INFINITY } else { 0.0 }
        } else {
            self.stats.transitions as f64 / self.stats.transversions as f64
        }
    }

    pub fn kappa_jukes_approx(&self) -> f64 {
        2.0 * self.ti_tv()
    }

    pub fn wobble_mod3_pct(&self) -> f64 {
        if self.stats.snps == 0 {
            0.0
        } else {
            (self.stats.wobble_mod3_hits as f64 * 100.0) / (self.stats.snps as f64)
        }
    }

    pub fn panel_conservation_pct(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            (self.stats.panel_match_bases as f64 * 100.0) / (self.total_bases as f64)
        }
    }

    pub fn avg_panel_tract_bp(&self) -> f64 {
        if self.stats.panel_tokens == 0 {
            0.0
        } else {
            self.stats.panel_match_bases as f64 / self.stats.panel_tokens as f64
        }
    }

    pub fn format_summary(&self) -> String {
        let mut s = format!(
            "=== EIMDL Compressed-Domain Evolutionary Report ===\n\
             Container     : {} (mode={}, units={})\n\
             Raw Bytes     : {}\n\
             Total Bases   : {} bp\n\
             Panel Copy    : {} bp ({:.3}%, {} tracts, avg {:.1} bp/tract)\n\
             Self / Tandem : {} bp / {} bp\n\
             Literals      : {} bp ({} runs, novelty={})\n\
             SNPs          : {} (Ti={}, Tv={}, Ti/Tv={:.4}, kappa~={:.4})\n\
             Codon Wobble  : {}/{} ({:.2}% mod-3 spaced)\n\
             Indel Splices : {} (INS={} bp, DEL={} bp)\n\
             Switches      : {} inter-donor recombination switches\n\
             Elapsed       : {:.4} s\n",
            self.format,
            self.mode,
            self.unit_count,
            self.raw_bytes,
            self.total_bases,
            self.stats.panel_match_bases,
            self.panel_conservation_pct(),
            self.stats.panel_tokens,
            self.avg_panel_tract_bp(),
            self.stats.self_match_bases,
            self.stats.tandem_match_bases,
            self.stats.literal_bases,
            self.stats.literal_runs,
            self.stats.novelty_len,
            self.stats.snps,
            self.stats.transitions,
            self.stats.transversions,
            self.ti_tv(),
            self.kappa_jukes_approx(),
            self.stats.wobble_mod3_hits,
            self.stats.snps,
            self.wobble_mod3_pct(),
            self.stats.splices,
            self.stats.inserted_bases,
            self.stats.deleted_bases,
            self.stats.recombination_switches,
            self.seconds,
        );
        if !self.donor_names.is_empty() {
            s.push_str("Donor Shares  :\n");
            let mut ranked: Vec<(usize, u64)> = self
                .donor_match_bases
                .iter()
                .enumerate()
                .filter(|&(_, &b)| b > 0)
                .map(|(i, &b)| (i, b))
                .collect();
            ranked.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
            for (id, bases) in ranked.into_iter().take(12) {
                let pct = if self.total_bases == 0 {
                    0.0
                } else {
                    (bases as f64 * 100.0) / (self.total_bases as f64)
                };
                s.push_str(&format!("  - {}: {} bp ({:.2}%)\n", self.donor_names[id], bases, pct));
            }
        }
        if let Some(ref c) = self.cds {
            s.push_str(&format!(
                "Proteomics    : {} ORFs (>=30aa), longest={} aa, total={} aa, GC={:.1}%\n",
                c.orf_count,
                c.longest_orf_aa,
                c.total_orf_aa,
                c.gc_permille as f64 * 0.1
            ));
        }
        s
    }

    pub fn to_json(&self) -> String {
        let mut out = String::from("{\n");
        out.push_str(&format!("  \"format\": \"{}\",\n", self.format));
        out.push_str(&format!("  \"is_fastq\": {},\n", self.is_fastq));
        out.push_str(&format!("  \"mode\": {},\n", self.mode));
        out.push_str(&format!("  \"raw_bytes\": {},\n", self.raw_bytes));
        out.push_str(&format!("  \"total_bases\": {},\n", self.total_bases));
        out.push_str(&format!("  \"panel_match_bases\": {},\n", self.stats.panel_match_bases));
        out.push_str(&format!("  \"self_match_bases\": {},\n", self.stats.self_match_bases));
        out.push_str(&format!("  \"tandem_match_bases\": {},\n", self.stats.tandem_match_bases));
        out.push_str(&format!("  \"literal_bases\": {},\n", self.stats.literal_bases));
        out.push_str(&format!("  \"panel_tokens\": {},\n", self.stats.panel_tokens));
        out.push_str(&format!("  \"snps\": {},\n", self.stats.snps));
        out.push_str(&format!("  \"transitions\": {},\n", self.stats.transitions));
        out.push_str(&format!("  \"transversions\": {},\n", self.stats.transversions));
        out.push_str(&format!("  \"ti_tv\": {:.6},\n", if self.ti_tv().is_finite() { self.ti_tv() } else { 0.0 }));
        out.push_str(&format!("  \"kappa_proxy\": {:.6},\n", if self.kappa_jukes_approx().is_finite() { self.kappa_jukes_approx() } else { 0.0 }));
        out.push_str(&format!("  \"wobble_mod3_hits\": {},\n", self.stats.wobble_mod3_hits));
        out.push_str(&format!("  \"wobble_mod3_pct\": {:.4},\n", self.wobble_mod3_pct()));
        out.push_str(&format!("  \"splices\": {},\n", self.stats.splices));
        out.push_str(&format!("  \"inserted_bases\": {},\n", self.stats.inserted_bases));
        out.push_str(&format!("  \"deleted_bases\": {},\n", self.stats.deleted_bases));
        out.push_str(&format!("  \"recombination_switches\": {},\n", self.stats.recombination_switches));
        out.push_str(&format!("  \"elapsed_seconds\": {:.6}\n", self.seconds));
        out.push('}');
        out
    }
}

struct LoadedContainer<'a> {
    format_name: &'static str,
    mode: u8,
    is_fastq: bool,
    raw_bytes: u64,
    panel_names: Vec<String>,
    panel_seqs: Vec<Vec<u8>>,
    nov_super: Vec<u8>,
    rank_to_unit: Vec<usize>,
    max_anc_rank: Vec<usize>,
    mosaic_slices: Vec<&'a [u8]>,
    cds: Option<CdsSummary>,
}

fn load_container<'a>(
    bytes: &'a [u8],
    panel: Option<&Path>,
    threads: i32,
    bounded: bool,
) -> Result<LoadedContainer<'a>, String> {
    if bytes.len() < 5 {
        return Err("Archive truncated".into());
    }
    if &bytes[..4] == MAGIC3 && (bytes[4] & hgs2::FLAG3_FASTQ) != 0 {
        let p = parse_hgs3_fastq(bytes, true)?;
        let panel_seqs = if p.mode == MODE_REF && !p.panel_names.is_empty() {
            verify_and_load_panel_for_decode(panel, &p.panel_names, &p.panel_sha)?
        } else {
            Vec::new()
        };
        let nov_super = if p.raw_bytes == 0 || p.total_records == 0 {
            Vec::new()
        } else {
            ffi::ans_novelty_decode(p.ans_nov_comp, p.nov_len, threads)?
        };
        let cds = if !p.cds_comp.is_empty() {
            Some(deserialize_cds_summary(p.cds_comp).map_err(|e| {
                format!("deserialize_cds_summary(len={}) failed: {e} head={:02X?}", p.cds_comp.len(), &p.cds_comp[..p.cds_comp.len().min(12)])
            })?)
        } else {
            None
        };
        let mosaic_slices = if p.raw_bytes == 0 || p.total_records == 0 {
            Vec::new()
        } else {
            vec![p.mosaic_bytes]
        };
        return Ok(LoadedContainer {
            format_name: "HGS3-FASTQ",
            mode: p.mode,
            is_fastq: true,
            raw_bytes: p.raw_bytes,
            panel_names: p.panel_names,
            panel_seqs,
            nov_super,
            rank_to_unit: vec![0; mosaic_slices.len()],
            max_anc_rank: vec![0; mosaic_slices.len()],
            mosaic_slices,
            cds,
        });
    }

    let expected_magic = if &bytes[..4] == MAGIC3 { MAGIC3 } else { MAGIC };
    let p: ParsedHgs2<'a> = parse_container(bytes, expected_magic, true)?;
    let panel_seqs = if p.mode == MODE_REF && !p.panel_names.is_empty() {
        verify_and_load_panel_for_decode(panel, &p.panel_names, &p.panel_sha)?
    } else {
        Vec::new()
    };
    let nov_super = if p.raw_len == 0 || p.num_units == 0 {
        Vec::new()
    } else if p.is_hgs3 {
        ffi::ans_novelty_decode(p.novelty_comp, p.total_novelty_bases, threads)?
    } else {
        ffi::novelty_decode(p.novelty_comp, p.total_novelty_bases, threads, bounded)?
    };
    let cds = if !p.cds_comp.is_empty() {
        Some(deserialize_cds_summary(p.cds_comp)?)
    } else {
        None
    };

    Ok(LoadedContainer {
        format_name: if p.is_hgs3 { "HGS3-FASTA" } else { "HGS2" },
        mode: p.mode,
        is_fastq: false,
        raw_bytes: p.raw_len,
        panel_names: p.panel_names,
        panel_seqs,
        nov_super,
        rank_to_unit: p.rank_to_unit,
        max_anc_rank: p.max_anc_rank,
        mosaic_slices: p.mosaic_slices,
        cds,
    })
}

type RankInspectBatch = (Vec<UnitInspectData>, Vec<(usize, usize)>);

fn inspect_all_ranks(
    c: &LoadedContainer<'_>,
    collect_variants: bool,
    collect_spans: bool,
) -> Result<RankInspectBatch, String> {
    let n = c.mosaic_slices.len();
    let mut reports = Vec::with_capacity(n);
    let mut nov_ranges = Vec::with_capacity(n);
    let mut nov_off = 0usize;

    if c.mode == MODE_COLLECTION && n > 1 {
        let needed_anc = c.max_anc_rank.iter().copied().max().unwrap_or(0);
        let mut anc_seqs: Vec<Vec<u8>> = Vec::with_capacity(needed_anc);
        for r in 0..n {
            let (_tl, nl) = ffi::mosaic_peek_lengths(c.mosaic_slices[r])?;
            let nstart = nov_off;
            let nend = nstart + (nl as usize);
            nov_off = nend;
            let nslice = &c.nov_super[nstart..nend];
            let anc = c.max_anc_rank[r];
            let donor_refs: Vec<&[u8]> = anc_seqs[..anc].iter().map(|v| v.as_slice()).collect();
            let rep = ffi::mosaic_inspect(
                c.mosaic_slices[r],
                nslice,
                &donor_refs,
                donor_refs.len(),
                collect_variants,
                collect_spans,
            )?;
            if r < needed_anc {
                let seq = ffi::decode_mosaic(c.mosaic_slices[r], nslice, &donor_refs)?;
                anc_seqs.push(seq);
            }
            reports.push(UnitInspectData {
                unit_index: c.rank_to_unit[r],
                stats: rep.stats,
                variants: rep.variants,
                spans: rep.spans,
                var_seq_pool: rep.var_seq_pool,
            });
            nov_ranges.push((nstart, nend));
        }
    } else {
        let panel_refs: Vec<&[u8]> = c.panel_seqs.iter().map(|v| v.as_slice()).collect();
        for u in 0..n {
            let (_tl, nl) = ffi::mosaic_peek_lengths(c.mosaic_slices[u])?;
            let nstart = nov_off;
            let nend = nstart + (nl as usize);
            nov_off = nend;
            let nslice = &c.nov_super[nstart..nend];
            let rep = ffi::mosaic_inspect(
                c.mosaic_slices[u],
                nslice,
                &panel_refs,
                panel_refs.len(),
                collect_variants,
                collect_spans,
            )?;
            reports.push(UnitInspectData {
                unit_index: u,
                stats: rep.stats,
                variants: rep.variants,
                spans: rep.spans,
                var_seq_pool: rep.var_seq_pool,
            });
            nov_ranges.push((nstart, nend));
        }
    }
    Ok((reports, nov_ranges))
}

fn accum_stats(dst: &mut HgsMosaicStats, src: &HgsMosaicStats) {
    dst.target_len += src.target_len;
    dst.novelty_len += src.novelty_len;
    dst.panel_match_bases += src.panel_match_bases;
    dst.self_match_bases += src.self_match_bases;
    dst.tandem_match_bases += src.tandem_match_bases;
    dst.literal_bases += src.literal_bases;
    dst.panel_tokens += src.panel_tokens;
    dst.self_tokens += src.self_tokens;
    dst.tandem_tokens += src.tandem_tokens;
    dst.literal_runs += src.literal_runs;
    dst.snps += src.snps;
    dst.transitions += src.transitions;
    dst.transversions += src.transversions;
    dst.wobble_mod3_hits += src.wobble_mod3_hits;
    dst.splices += src.splices;
    dst.inserted_bases += src.inserted_bases;
    dst.deleted_bases += src.deleted_bases;
    dst.recombination_switches += src.recombination_switches;
}

pub fn analyze_archive(
    archive: &Path,
    panel: Option<&Path>,
    threads: i32,
    bounded: bool,
) -> Result<AnalyzeReport, String> {
    let t0 = Instant::now();
    let bytes = std::fs::read(archive).map_err(|e| e.to_string())?;
    let c = load_container(&bytes, panel, threads, bounded)?;
    let (units, _) = inspect_all_ranks(&c, false, !c.panel_names.is_empty())?;

    let mut total_stats = HgsMosaicStats::default();
    let mut donor_bases = vec![0u64; c.panel_names.len()];
    for u in &units {
        accum_stats(&mut total_stats, &u.stats);
        for sp in &u.spans {
            if sp.kind == 1 && (sp.donor_id as usize) < donor_bases.len() {
                donor_bases[sp.donor_id as usize] += sp.tgt_len;
            }
        }
    }

    Ok(AnalyzeReport {
        format: c.format_name,
        is_fastq: c.is_fastq,
        mode: c.mode,
        raw_bytes: c.raw_bytes,
        total_bases: total_stats.target_len,
        stats: total_stats,
        donor_names: c.panel_names,
        donor_match_bases: donor_bases,
        unit_count: units.len(),
        cds: c.cds,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

pub struct ExtractedVcf {
    pub vcf_text: String,
    pub total_variants: usize,
    pub snvs: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub complex_indels: usize,
    pub literal_runs_skipped: usize,
    pub read_boundary_filtered: usize,
    pub seconds: f64,
}

fn append_two_bit_ascii(out: &mut String, seq: &[u8]) {
    const MAP: [char; 4] = ['A', 'C', 'G', 'T'];
    for &b in seq {
        out.push(MAP[(b & 3) as usize]);
    }
}

fn crosses_fastq_read_edge(tgt_pos: u64, alt_len: u32, fixed_len: usize, read_ends: &[u64]) -> bool {
    let span = u64::from(alt_len.max(1));
    if fixed_len > 0 {
        let r = fixed_len as u64;
        let r_start = (tgt_pos / r) * r;
        let r_end = r_start + r;
        tgt_pos < r_start + 8 || tgt_pos.saturating_add(span).saturating_add(8) > r_end
    } else if !read_ends.is_empty() {
        let idx = read_ends.partition_point(|&e| e <= tgt_pos);
        let r_start = if idx == 0 { 0 } else { read_ends[idx - 1] };
        let r_end = *read_ends.get(idx).unwrap_or(&u64::MAX);
        tgt_pos < r_start + 8 || tgt_pos.saturating_add(span).saturating_add(8) > r_end
    } else {
        false
    }
}

pub fn extract_variants_vcf(
    archive: &Path,
    panel: Option<&Path>,
    out_vcf: Option<&Path>,
    threads: i32,
    bounded: bool,
) -> Result<ExtractedVcf, String> {
    extract_variants_vcf_with_support(archive, panel, out_vcf, threads, bounded, 1)
}

pub fn extract_variants_vcf_with_support(
    archive: &Path,
    panel: Option<&Path>,
    out_vcf: Option<&Path>,
    threads: i32,
    bounded: bool,
    min_support: usize,
) -> Result<ExtractedVcf, String> {
    use std::collections::BTreeMap;

    let t0 = Instant::now();
    let bytes = std::fs::read(archive).map_err(|e| e.to_string())?;
    let c = load_container(&bytes, panel, threads, bounded)?;
    let (units, _) = inspect_all_ranks(&c, true, false)?;

    let (fq_fixed, fq_ends) = if c.is_fastq && bytes.len() >= 5 && &bytes[..4] == MAGIC3 {
        let pf = parse_hgs3_fastq(&bytes, false)?;
        hgs3::extract_fastq_read_bounds(pf.meta_comp, pf.total_records as usize)
            .unwrap_or((pf.read_len_hint as usize, Vec::new()))
    } else {
        (0usize, Vec::new())
    };

    struct AggCall {
        kind: u8,
        kind_str: &'static str,
        tgt_pos: u64,
        unit_index: usize,
        strand: u8,
        is_transition: u8,
        prev_mod3: u8,
        ac: usize,
    }

    let mut agg: BTreeMap<(usize, u64, Vec<u8>, Vec<u8>), AggCall> = BTreeMap::new();
    let mut literal_runs_skipped = 0usize;
    let mut read_boundary_filtered = 0usize;

    for u in &units {
        for v in &u.variants {
            if v.kind == 4 {
                literal_runs_skipped += 1;
                continue;
            }
            if c.is_fastq && v.kind != 0 && crosses_fastq_read_edge(v.tgt_pos, v.alt_len, fq_fixed, &fq_ends) {
                read_boundary_filtered += 1;
                continue;
            }
            let did = v.donor_id as usize;
            if did >= c.panel_seqs.len() {
                continue;
            }
            let donor = &c.panel_seqs[did];
            let dlen = donor.len() as u64;

            let (fwd_start, fwd_end) = if v.strand == 0 {
                (v.ref_pos, v.ref_pos + u64::from(v.ref_len))
            } else {
                let end = dlen.saturating_sub(v.ref_pos);
                let start = end.saturating_sub(u64::from(v.ref_len));
                (start, end)
            };
            if fwd_end > dlen {
                continue;
            }

            let mut fwd_ref = donor[fwd_start as usize..fwd_end as usize].to_vec();
            let mut fwd_alt = if v.alt_len == 0 || v.alt_seq_off == u64::MAX {
                Vec::new()
            } else {
                let off = v.alt_seq_off as usize;
                let slice = &u.var_seq_pool[off..off + (v.alt_len as usize)];
                if v.strand == 0 {
                    slice.to_vec()
                } else {
                    slice.iter().rev().map(|&b| (b ^ 3) & 3).collect()
                }
            };

            let pos1: u64 = if fwd_ref.is_empty() || fwd_alt.is_empty() {
                if fwd_start > 0 {
                    let anc = donor[(fwd_start - 1) as usize] & 3;
                    fwd_ref.insert(0, anc);
                    fwd_alt.insert(0, anc);
                    fwd_start
                } else if fwd_end < dlen {
                    let anc = donor[fwd_end as usize] & 3;
                    fwd_ref.push(anc);
                    fwd_alt.push(anc);
                    1
                } else {
                    fwd_start + 1
                }
            } else {
                fwd_start + 1
            };

            let kind_str = match v.kind {
                0 => "SNV",
                1 => "INS",
                2 => "DEL",
                _ => "COMPLEX",
            };

            agg.entry((did, pos1, fwd_ref, fwd_alt))
                .and_modify(|e| e.ac += 1)
                .or_insert(AggCall {
                    kind: v.kind,
                    kind_str,
                    tgt_pos: v.tgt_pos,
                    unit_index: u.unit_index,
                    strand: v.strand,
                    is_transition: v.is_transition,
                    prev_mod3: v.prev_mod3,
                    ac: 1,
                });
        }
    }

    let mut vcf = String::with_capacity(4096);
    vcf.push_str("##fileformat=VCFv4.2\n");
    vcf.push_str("##source=EIMDL-HyperGraSS-CompressedDomain\n");
    vcf.push_str("##INFO=<ID=TGT_POS,Number=1,Type=Integer,Description=\"0-based coordinate in target archive\">\n");
    vcf.push_str("##INFO=<ID=UNIT,Number=1,Type=Integer,Description=\"Target archive unit index\">\n");
    vcf.push_str("##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Donor alignment strand (+ or -)\">\n");
    vcf.push_str("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant class (SNV, INS, DEL, COMPLEX)\">\n");
    vcf.push_str("##INFO=<ID=TI,Number=1,Type=Integer,Description=\"1 if transition SNP else 0\">\n");
    vcf.push_str("##INFO=<ID=MOD3,Number=1,Type=Integer,Description=\"1 if SNP step % 3 == 0 codon-wobble spaced\">\n");
    vcf.push_str("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Supporting allele observations across units/reads\">\n");
    for (i, dname) in c.panel_names.iter().enumerate() {
        let dlen = c.panel_seqs.get(i).map_or(0, |v| v.len());
        vcf.push_str(&format!("##contig=<ID={dname},length={dlen}>\n"));
    }
    vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    let min_sup = min_support.max(1);
    let mut total_variants = 0usize;
    let mut snvs = 0usize;
    let mut insertions = 0usize;
    let mut deletions = 0usize;
    let mut complex_indels = 0usize;

    for ((did, pos1, fwd_ref, fwd_alt), call) in agg {
        if call.ac < min_sup {
            continue;
        }
        match call.kind {
            0 => snvs += 1,
            1 => insertions += 1,
            2 => deletions += 1,
            _ => complex_indels += 1,
        }
        total_variants += 1;
        let chrom = &c.panel_names[did];
        let mut ref_str = String::with_capacity(fwd_ref.len());
        let mut alt_str = String::with_capacity(fwd_alt.len());
        append_two_bit_ascii(&mut ref_str, &fwd_ref);
        append_two_bit_ascii(&mut alt_str, &fwd_alt);
        let strand_ch = if call.strand == 0 { '+' } else { '-' };

        vcf.push_str(&format!(
            "{chrom}\t{pos1}\t.\t{ref_str}\t{alt_str}\t.\tPASS\tTGT_POS={};UNIT={};STRAND={strand_ch};TYPE={};TI={};MOD3={};AC={}\n",
            call.tgt_pos, call.unit_index, call.kind_str, call.is_transition, call.prev_mod3, call.ac
        ));
    }

    if let Some(p) = out_vcf {
        let mut w = BufWriter::new(File::create(p).map_err(|e| e.to_string())?);
        w.write_all(vcf.as_bytes()).map_err(|e| e.to_string())?;
        w.flush().map_err(|e| e.to_string())?;
    }

    Ok(ExtractedVcf {
        vcf_text: vcf,
        total_variants,
        snvs,
        insertions,
        deletions,
        complex_indels,
        literal_runs_skipped,
        read_boundary_filtered,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

fn parse_ascii_query(q: &str) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(q.len());
    for &b in q.as_bytes() {
        match b {
            b'A' | b'a' => out.push(0),
            b'C' | b'c' => out.push(1),
            b'G' | b'g' => out.push(2),
            b'T' | b't' | b'U' | b'u' => out.push(3),
            _ => return Err(format!("Unsupported non-ACGT base '{}' in query", b as char)),
        }
    }
    if out.is_empty() {
        return Err("Query sequence must be non-empty".into());
    }
    Ok(out)
}

fn revcomp_2bit(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| (b ^ 3) & 3).collect()
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct KmerHit {
    pub unit_index: usize,
    pub tgt_pos: u64,
    pub strand: u8,
}

pub struct KmerSearchReport {
    pub query_len: usize,
    pub hits: Vec<KmerHit>,
    pub exact_zero_target_spans: usize,
    pub patched_or_local_spans: usize,
    pub seconds: f64,
}

#[inline]
fn donor_base_strand(donor: &[u8], strand: u8, mr: usize) -> u8 {
    if strand == 0 {
        donor[mr] & 3
    } else {
        (donor[donor.len() - 1 - mr] ^ 3) & 3
    }
}

fn resolve_unit_base_from_spans(
    pos: usize,
    spans: &[HgsSpanEvent],
    snps: &[&HgsVariantRecord],
    var_seq_pool: &[u8],
    novelty: &[u8],
    donors: &[&[u8]],
) -> u8 {
    let mut cur_pos = pos;
    let mut flip = 0u8;
    for _depth in 0..64 {
        let sidx = spans
            .partition_point(|s| (s.tgt_pos as usize) <= cur_pos)
            .saturating_sub(1);
        let sp = &spans[sidx];
        let off = cur_pos - (sp.tgt_pos as usize);
        match sp.kind {
            0 => return (novelty[sp.ref_pos as usize + off] & 3) ^ flip,
            1 => {
                if sp.has_edits != 0
                    && let Ok(vidx) = snps.binary_search_by_key(&(cur_pos as u64), |v| v.tgt_pos)
                {
                    let vr = snps[vidx];
                    if vr.alt_seq_off != u64::MAX {
                        return (var_seq_pool[vr.alt_seq_off as usize] & 3) ^ flip;
                    }
                }
                return donor_base_strand(donors[sp.donor_id as usize], sp.strand, sp.ref_pos as usize + off) ^ flip;
            }
            3 => {
                let period = (sp.ref_pos as usize).max(1);
                cur_pos = (sp.tgt_pos as usize) - period + (off % period);
            }
            _ => {
                let mr = sp.ref_pos as usize + off;
                if sp.strand == 0 {
                    cur_pos = mr;
                } else {
                    cur_pos = (sp.tgt_pos as usize) - 1 - mr;
                    flip ^= 3;
                }
            }
        }
    }
    0
}

fn resolve_span_pos(
    pos: usize,
    sp: &HgsSpanEvent,
    u: &UnitInspectData,
    snps: &[&HgsVariantRecord],
    novelty: &[u8],
    donors: &[&[u8]],
) -> u8 {
    let off = pos - (sp.tgt_pos as usize);
    match sp.kind {
        0 => novelty[sp.ref_pos as usize + off] & 3,
        1 => {
            if sp.has_edits != 0
                && let Ok(vi) = snps.binary_search_by_key(&(pos as u64), |v| v.tgt_pos)
            {
                let vr = snps[vi];
                if vr.alt_seq_off != u64::MAX {
                    return u.var_seq_pool[vr.alt_seq_off as usize] & 3;
                }
            }
            donor_base_strand(donors[sp.donor_id as usize], sp.strand, sp.ref_pos as usize + off)
        }
        _ => resolve_unit_base_from_spans(pos, &u.spans, snps, &u.var_seq_pool, novelty, donors),
    }
}

fn scan_window_both_strands(
    win: &[u8],
    base_tgt_pos: u64,
    unit_idx: usize,
    q_fwd: &[u8],
    q_rc: Option<&[u8]>,
    hits: &mut Vec<KmerHit>,
) {
    let k = q_fwd.len();
    if win.len() < k {
        return;
    }
    for i in 0..=win.len() - k {
        if &win[i..i + k] == q_fwd {
            hits.push(KmerHit {
                unit_index: unit_idx,
                tgt_pos: base_tgt_pos + i as u64,
                strand: 0,
            });
        }
        if let Some(rc) = q_rc
            && &win[i..i + k] == rc
        {
            hits.push(KmerHit {
                unit_index: unit_idx,
                tgt_pos: base_tgt_pos + i as u64,
                strand: 1,
            });
        }
    }
}

pub fn search_kmers(
    archive: &Path,
    panel: Option<&Path>,
    query_ascii: &str,
    both_strands: bool,
    threads: i32,
    bounded: bool,
) -> Result<KmerSearchReport, String> {
    let t0 = Instant::now();
    let q_fwd = parse_ascii_query(query_ascii)?;
    let q_rc = if both_strands {
        Some(revcomp_2bit(&q_fwd))
    } else {
        None
    };
    let k = q_fwd.len();

    let bytes = std::fs::read(archive).map_err(|e| e.to_string())?;
    let c = load_container(&bytes, panel, threads, bounded)?;
    if c.mode == MODE_COLLECTION && c.mosaic_slices.len() > 1 {
        // Materialize unit ranks cleanly in Mode C DAG for exact donor lookup
        let n = c.mosaic_slices.len();
        let mut nov_off = 0usize;
        let mut decoded_ranks: Vec<Vec<u8>> = Vec::with_capacity(n);
        let mut hits = Vec::new();
        let mut patched = 0usize;
        for r in 0..n {
            let (_tl, nl) = ffi::mosaic_peek_lengths(c.mosaic_slices[r])?;
            let nslice = &c.nov_super[nov_off..nov_off + nl as usize];
            nov_off += nl as usize;
            let anc = c.max_anc_rank[r];
            let donor_refs: Vec<&[u8]> = decoded_ranks[..anc].iter().map(|v| v.as_slice()).collect();
            let dna = ffi::decode_mosaic(c.mosaic_slices[r], nslice, &donor_refs)?;
            scan_window_both_strands(&dna, 0, c.rank_to_unit[r], &q_fwd, q_rc.as_deref(), &mut hits);
            patched += 1;
            decoded_ranks.push(dna);
        }
        hits.sort();
        return Ok(KmerSearchReport {
            query_len: k,
            hits,
            exact_zero_target_spans: 0,
            patched_or_local_spans: patched,
            seconds: t0.elapsed().as_secs_f64(),
        });
    }

    let panel_refs: Vec<&[u8]> = c.panel_seqs.iter().map(|v| v.as_slice()).collect();
    let (units, nov_ranges) = inspect_all_ranks(&c, true, true)?;
    let mut hits = Vec::new();
    let mut exact_zero_target_spans = 0usize;
    let mut patched_or_local_spans = 0usize;

    for (idx, u) in units.iter().enumerate() {
        let (ns, ne) = nov_ranges[idx];
        let novelty = &c.nov_super[ns..ne];
        let has_self_edits = u.spans.iter().any(|s| s.kind >= 2 && s.has_edits != 0);
        let mut lazy_self_dna: Option<Vec<u8>> = None;

        let mut snps: Vec<&HgsVariantRecord> = u.variants.iter().filter(|v| v.kind == 0).collect();
        snps.sort_by_key(|v| v.tgt_pos);
        // Resolve SNP alt base directly into a lightweight lookup helper
        let mut snp_cursor = 0usize;
        let mut span_buf: Vec<u8> = Vec::new();
        let mut junction_buf: Vec<u8> = Vec::new();
        let mut prev_tail: Vec<u8> = Vec::new();

        for sp in &u.spans {
            let sp_start = sp.tgt_pos as usize;
            let sp_len = sp.tgt_len as usize;
            let sp_end = sp_start + sp_len;

            if sp.kind >= 2 && has_self_edits && lazy_self_dna.is_none() {
                lazy_self_dna = Some(ffi::decode_mosaic(c.mosaic_slices[idx], novelty, &panel_refs)?);
            }

            // Advance SNP cursor past prior spans
            while snp_cursor < snps.len() && (snps[snp_cursor].tgt_pos as usize) < sp_start {
                snp_cursor += 1;
            }

            if !prev_tail.is_empty() && k > 1 {
                let take_right = sp_len.min(k - 1);
                junction_buf.clear();
                junction_buf.extend_from_slice(&prev_tail);
                for off in 0..take_right {
                    let b = if sp.kind >= 2 && let Some(ref sdna) = lazy_self_dna {
                        sdna[sp_start + off] & 3
                    } else {
                        resolve_span_pos(sp_start + off, sp, u, &snps, novelty, &panel_refs)
                    };
                    junction_buf.push(b);
                }
                let jstart = sp_start - prev_tail.len();
                let limit_start = prev_tail.len();
                if junction_buf.len() >= k {
                    for i in 0..limit_start.min(junction_buf.len() - k + 1) {
                        if &junction_buf[i..i + k] == q_fwd.as_slice() {
                            hits.push(KmerHit {
                                unit_index: u.unit_index,
                                tgt_pos: (jstart + i) as u64,
                                strand: 0,
                            });
                        }
                        if let Some(ref rc) = q_rc
                            && &junction_buf[i..i + k] == rc.as_slice()
                        {
                            hits.push(KmerHit {
                                unit_index: u.unit_index,
                                tgt_pos: (jstart + i) as u64,
                                strand: 1,
                            });
                        }
                    }
                }
            }

            // Span-interior exact search
            if sp_len >= k {
                if sp.kind == 0 {
                    exact_zero_target_spans += 1;
                    let lit = &novelty[sp.ref_pos as usize..sp.ref_pos as usize + sp_len];
                    scan_window_both_strands(lit, sp.tgt_pos, u.unit_index, &q_fwd, q_rc.as_deref(), &mut hits);
                } else if sp.kind == 1 && sp.has_edits == 0 && sp.strand == 0 {
                    exact_zero_target_spans += 1;
                    let d = panel_refs[sp.donor_id as usize];
                    let rs = sp.ref_pos as usize;
                    scan_window_both_strands(&d[rs..rs + sp_len], sp.tgt_pos, u.unit_index, &q_fwd, q_rc.as_deref(), &mut hits);
                } else if sp.kind == 1 && sp.has_edits == 0 && sp.strand == 1 {
                    exact_zero_target_spans += 1;
                    let d = panel_refs[sp.donor_id as usize];
                    let dlen = d.len();
                    let rs = sp.ref_pos as usize;
                    for off in 0..=sp_len - k {
                        let mut mf = true;
                        let mut mr = q_rc.is_some();
                        for j in 0..k {
                            let db = (d[dlen - 1 - (rs + off + j)] ^ 3) & 3;
                            if db != q_fwd[j] {
                                mf = false;
                            }
                            if let Some(ref rc) = q_rc
                                && db != rc[j]
                            {
                                mr = false;
                            }
                            if !mf && !mr {
                                break;
                            }
                        }
                        if mf {
                            hits.push(KmerHit { unit_index: u.unit_index, tgt_pos: sp.tgt_pos + off as u64, strand: 0 });
                        }
                        if mr {
                            hits.push(KmerHit { unit_index: u.unit_index, tgt_pos: sp.tgt_pos + off as u64, strand: 1 });
                        }
                    }
                } else if sp.kind >= 2 && let Some(ref sdna) = lazy_self_dna {
                    patched_or_local_spans += 1;
                    scan_window_both_strands(&sdna[sp_start..sp_end], sp.tgt_pos, u.unit_index, &q_fwd, q_rc.as_deref(), &mut hits);
                } else {
                    patched_or_local_spans += 1;
                    span_buf.resize(sp_len, 0);
                    if sp.kind == 1 {
                        let d = panel_refs[sp.donor_id as usize];
                        let rs = sp.ref_pos as usize;
                        for (off, dst) in span_buf.iter_mut().enumerate() {
                            *dst = donor_base_strand(d, sp.strand, rs + off);
                        }
                        let mut sc = snp_cursor;
                        while sc < snps.len() && (snps[sc].tgt_pos as usize) < sp_end {
                            let pos = snps[sc].tgt_pos as usize;
                            span_buf[pos - sp_start] = u.var_seq_pool[snps[sc].alt_seq_off as usize] & 3;
                            sc += 1;
                        }
                    } else {
                        for (off, dst) in span_buf.iter_mut().enumerate() {
                            *dst = resolve_unit_base_from_spans(
                                sp_start + off,
                                &u.spans,
                                &snps,
                                &u.var_seq_pool,
                                novelty,
                                &panel_refs,
                            );
                        }
                    }
                    scan_window_both_strands(&span_buf, sp.tgt_pos, u.unit_index, &q_fwd, q_rc.as_deref(), &mut hits);
                }
            } else {
                exact_zero_target_spans += 1;
            }

            // Maintain trailing (k - 1) target bases cleanly without buffering whole target
            if k > 1 {
                let keep_prev = (k - 1).saturating_sub(sp_len);
                let tail_start = if keep_prev > 0 { prev_tail.len().saturating_sub(keep_prev) } else { prev_tail.len() };
                let mut next_tail = prev_tail[tail_start..].to_vec();
                let take_this = sp_len.min((k - 1) - next_tail.len());
                let start_off = sp_len - take_this;
                for off in start_off..sp_len {
                    let b = if sp.kind >= 2 && let Some(ref sdna) = lazy_self_dna {
                        sdna[sp_start + off] & 3
                    } else {
                        resolve_span_pos(sp_start + off, sp, u, &snps, novelty, &panel_refs)
                    };
                    next_tail.push(b);
                }
                prev_tail = next_tail;
            }
        }
    }

    hits.sort();
    hits.dedup();

    Ok(KmerSearchReport {
        query_len: k,
        hits,
        exact_zero_target_spans,
        patched_or_local_spans,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

pub fn query_slice(
    archive: &Path,
    panel: Option<&Path>,
    unit_idx: usize,
    start: u64,
    end: u64,
    threads: i32,
    bounded: bool,
) -> Result<String, String> {
    let bytes = std::fs::read(archive).map_err(|e| e.to_string())?;
    let c = load_container(&bytes, panel, threads, bounded)?;
    let n = c.mosaic_slices.len();
    if unit_idx >= n {
        return Err(format!("Unit index {unit_idx} out of range 0..{n}"));
    }

    let two_bit = if c.mode == MODE_COLLECTION && n > 1 {
        let mut nov_off = 0usize;
        let mut decoded: Vec<Vec<u8>> = Vec::with_capacity(n);
        let mut found = Vec::new();
        for r in 0..n {
            let (_tl, nl) = ffi::mosaic_peek_lengths(c.mosaic_slices[r])?;
            let nslice = &c.nov_super[nov_off..nov_off + nl as usize];
            nov_off += nl as usize;
            let anc = c.max_anc_rank[r];
            let drefs: Vec<&[u8]> = decoded[..anc].iter().map(|v| v.as_slice()).collect();
            let seq = ffi::decode_mosaic(c.mosaic_slices[r], nslice, &drefs)?;
            if c.rank_to_unit[r] == unit_idx {
                found = seq;
                break;
            }
            decoded.push(seq);
        }
        found
    } else {
        let mut nov_off = 0usize;
        for mos in c.mosaic_slices.iter().take(unit_idx) {
            let (_tl, nl) = ffi::mosaic_peek_lengths(mos)?;
            nov_off += nl as usize;
        }
        let (_tl, nl) = ffi::mosaic_peek_lengths(c.mosaic_slices[unit_idx])?;
        let nslice = &c.nov_super[nov_off..nov_off + nl as usize];
        let panel_refs: Vec<&[u8]> = c.panel_seqs.iter().map(|v| v.as_slice()).collect();
        ffi::decode_mosaic(c.mosaic_slices[unit_idx], nslice, &panel_refs)?
    };

    let s = (start as usize).min(two_bit.len());
    let e = (end as usize).max(s).min(two_bit.len());
    let mut out = String::with_capacity(e - s);
    append_two_bit_ascii(&mut out, &two_bit[s..e]);
    Ok(out)
}
