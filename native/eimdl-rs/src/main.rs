use std::env;
use std::fs;
use std::path::PathBuf;
use std::process;

use eimdl_native::compute;
use eimdl_native::ffi;
use eimdl_native::hgs2::{self, CompressOptions, MODE_COLLECTION, MODE_SINGLE};
use eimdl_native::hgs3::{self, Hgs3CompressOptions, OmicInputKind};

fn print_usage() {
    eprintln!(
        "EIMDL / HyperGraSS v3.0.0 Native Engine (C23 AVX2/4x-rANS + C++23 64-bit DAG/CRM/Inspect + Rust 2024)\n\
         Usage:\n\
           eimdl --version\n\
           eimdl compress --target <file> --out <archive.hgs2> [--panel <dir|file>] [--collection] [--single] [--bounded] [--threads <N>]\n\
           eimdl decompress --archive <archive> [--panel <dir|file>] [--bounded] [--threads <N>] --out <file>\n\
           eimdl compress-omics --target <fasta|fastq> --out <archive.hgs3> [--format auto|fasta|fastq] [--panel <dir|file>] [--cds] [--threads <N>]\n\
           eimdl decompress-omics --archive <archive.hgs3|hgs2> [--panel <dir|file>] [--threads <N>] --out <file>\n\
           eimdl analyze --archive <archive.hgs2|hgs3> [--panel <dir|file>] [--json] [--threads <N>]\n\
           eimdl variants --archive <archive.hgs2|hgs3> --panel <dir|file> [--out <out.vcf>] [--json]\n\
           eimdl kmer --archive <archive.hgs2|hgs3> [--panel <dir|file>] --query <ACGT> [--both-strands] [--json]\n\
           eimdl query --archive <archive.hgs2|hgs3> [--panel <dir|file>] [--unit <N>] --start <pos> --end <pos>\n\
           eimdl inspect --archive <archive.hgs2>"
    );
}

fn dispatch_subcommand(sub: &str, rest: &[String]) -> Result<(), String> {
    match sub {
        "--version" | "-V" | "version" => {
            println!("hypergrass-eimdl-native 3.0.0 (avx2={})", ffi::has_avx2());
            Ok(())
        }
        "compress" => run_compress(rest),
        "decompress" | "decompress-omics" => run_decompress(rest),
        "compress-omics" => run_compress_omics(rest),
        "analyze" => run_analyze(rest),
        "variants" => run_variants(rest),
        "kmer" => run_kmer(rest),
        "query" => run_query(rest),
        "inspect" => run_inspect(rest),
        _ => {
            print_usage();
            Err(format!("Unknown subcommand: {sub}"))
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        print_usage();
        process::exit(1);
    }
    if let Err(e) = dispatch_subcommand(&args[1], &args[2..]) {
        eprintln!("Error: {e}");
        process::exit(1);
    }
}

fn next_val(args: &[String], i: &mut usize, flag: &str) -> Result<String, String> {
    *i += 1;
    if *i >= args.len() {
        return Err(format!("{flag} requires an argument"));
    }
    Ok(args[*i].clone())
}

fn parse_common_compress_flag(
    args: &[String],
    i: &mut usize,
    opts: &mut CompressOptions,
    target: &mut Option<PathBuf>,
    out: &mut Option<PathBuf>,
) -> Result<bool, String> {
    match args[*i].as_str() {
        "--target" | "--input" => *target = Some(PathBuf::from(next_val(args, i, "--target")?)),
        "--out" | "--output" => *out = Some(PathBuf::from(next_val(args, i, "--out")?)),
        "--panel" | "--ref" => opts.reference_panel = Some(PathBuf::from(next_val(args, i, "--panel")?)),
        "--exclude" => { opts.exclude_panel_names.insert(next_val(args, i, "--exclude")?); }
        "--panel-list" => {
            let lp = PathBuf::from(next_val(args, i, "--panel-list")?);
            let txt = fs::read_to_string(lp).map_err(|e| e.to_string())?;
            let list: Vec<String> = txt.lines().map(|l| l.trim().to_string()).filter(|l| !l.is_empty()).collect();
            opts.explicit_panel_names = Some(list);
        }
        "--panel-limit" => {
            opts.panel_limit = next_val(args, i, "--panel-limit")?.parse().map_err(|_| "Invalid --panel-limit".to_string())?;
        }
        "--threads" => {
            opts.threads = next_val(args, i, "--threads")?.parse::<i32>().map_err(|_| "Invalid --threads".to_string())?.max(1);
        }
        "--bounded" => opts.bounded = true,
        "--collection" => opts.mode = MODE_COLLECTION as i32,
        "--single" => { opts.single_genome = true; opts.mode = MODE_SINGLE as i32; }
        "--k" | "--step" => { let _ = next_val(args, i, "--k")?; }
        "--no-approx" => opts.mosaic_flags.allow_approx = false,
        "--no-splice" => opts.mosaic_flags.allow_splice = false,
        "--no-v4" => opts.mosaic_flags.allow_kimura = false,
        "--no-wobble" => opts.mosaic_flags.allow_wobble = false,
        "--no-self" => opts.mosaic_flags.allow_self = false,
        "--no-tandem" => opts.mosaic_flags.allow_tandem = false,
        _ => return Ok(false),
    }
    Ok(true)
}

fn run_compress(args: &[String]) -> Result<(), String> {
    let mut target: Option<PathBuf> = None;
    let mut out: Option<PathBuf> = None;
    let mut opts = CompressOptions::default();

    let mut i = 0;
    while i < args.len() {
        if !parse_common_compress_flag(args, &mut i, &mut opts, &mut target, &mut out)? {
            return Err(format!("Unknown compress option: {}", args[i]));
        }
        i += 1;
    }

    let tgt = target.ok_or("compress requires --target <path> and --out <path>")?;
    let out_p = out.ok_or("compress requires --target <path> and --out <path>")?;
    let rep = hgs2::compress_file(&tgt, &out_p, &opts)?;

    println!(
        "{{\"mode\":\"{}\",\"engine\":\"native-c-cpp-rust\",\"compress_seconds\":{:.6},\"raw_bytes\":{},\"target_bases\":{},\"archive_bytes\":{},\"bpb\":{:.6},\"ratio\":{:.4},\"panel_genomes\":{},\"panel_bases\":{},\"panel_match_bases\":{},\"self_match_bases\":{},\"tandem_match_bases\":{},\"novelty_bases\":{},\"panel_tokens\":{},\"self_tokens\":{},\"tandem_tokens\":{},\"literal_runs\":{},\"indel_splices\":{},\"substitutions\":{}}}",
        rep.mode_name(),
        rep.seconds,
        rep.raw_bytes,
        rep.target_dna_bases,
        rep.archive_bytes,
        rep.bpb(),
        rep.ratio(),
        rep.panel_genomes,
        rep.panel_bases,
        rep.panel_match_bases,
        rep.self_match_bases,
        rep.tandem_match_bases,
        rep.novelty_bases,
        rep.panel_tokens,
        rep.self_tokens,
        rep.tandem_tokens,
        rep.literal_runs,
        rep.indel_splices,
        rep.substitutions,
    );
    Ok(())
}

fn run_compress_omics(args: &[String]) -> Result<(), String> {
    let mut target: Option<PathBuf> = None;
    let mut out: Option<PathBuf> = None;
    let mut opts = Hgs3CompressOptions::default();

    let mut i = 0;
    while i < args.len() {
        if parse_common_compress_flag(args, &mut i, &mut opts.base, &mut target, &mut out)? {
            i += 1;
            continue;
        }
        match args[i].as_str() {
            "--format" => {
                opts.format = match next_val(args, &mut i, "--format")?.to_ascii_lowercase().as_str() {
                    "auto" => OmicInputKind::Auto,
                    "fasta" | "fa" | "fna" => OmicInputKind::Fasta,
                    "fastq" | "fq" => OmicInputKind::Fastq,
                    other => return Err(format!("Unsupported --format {other}")),
                };
            }
            "--cds" => opts.emit_cds_digest = true,
            other => return Err(format!("Unknown compress-omics option: {other}")),
        }
        i += 1;
    }

    let tgt = target.ok_or("compress-omics requires --target <path> and --out <path>")?;
    let out_p = out.ok_or("compress-omics requires --target <path> and --out <path>")?;
    let rep = hgs3::compress_omics(&tgt, &out_p, &opts)?;
    println!(
        "{{\"container\":\"HGS3\",\"is_fastq\":{},\"mode\":{},\"engine\":\"native-c-cpp-rust-ans4\",\"compress_seconds\":{:.6},\"raw_bytes\":{},\"archive_bytes\":{},\"bpb\":{:.6},\"ratio\":{:.4},\"records\":{},\"dna_bases\":{},\"qual_bytes\":{},\"meta_comp_bytes\":{},\"mosaic_comp_bytes\":{},\"novelty_comp_bytes\":{},\"qual_comp_bytes\":{},\"cds_orfs\":{}}}",
        rep.is_fastq,
        rep.mode,
        rep.seconds,
        rep.raw_bytes,
        rep.archive_bytes,
        rep.bpb(),
        rep.ratio(),
        rep.total_records,
        rep.dna_bases,
        rep.qual_bytes,
        rep.meta_comp_bytes,
        rep.mosaic_comp_bytes,
        rep.novelty_comp_bytes,
        rep.qual_comp_bytes,
        rep.cds.as_ref().map_or(0, |c| c.orf_count),
    );
    Ok(())
}

fn run_decompress(args: &[String]) -> Result<(), String> {
    let mut archive: Option<PathBuf> = None;
    let mut out: Option<PathBuf> = None;
    let mut panel: Option<PathBuf> = None;
    let mut bounded = false;
    let mut threads = std::thread::available_parallelism().map_or(4, |n| n.get() as i32);

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--archive" | "--input" => archive = Some(PathBuf::from(next_val(args, &mut i, "--archive")?)),
            "--panel" | "--ref" => panel = Some(PathBuf::from(next_val(args, &mut i, "--panel")?)),
            "--out" | "--output" => out = Some(PathBuf::from(next_val(args, &mut i, "--out")?)),
            "--bounded" => bounded = true,
            "--threads" => { threads = next_val(args, &mut i, "--threads")?.parse::<i32>().map_err(|_| "Invalid --threads".to_string())?.max(1); }
            other => return Err(format!("Unknown decompress option: {other}")),
        }
        i += 1;
    }

    let arc = archive.ok_or("decompress requires --archive <path> and --out <path>")?;
    let out_p = out.ok_or("decompress requires --archive <path> and --out <path>")?;
    let secs = hgs3::decompress_omics(&arc, &out_p, panel.as_deref(), threads, bounded)?;
    println!("{{\"decompress_seconds\":{secs:.6},\"engine\":\"native-c-cpp-rust\",\"verified_crc32\":true}}");
    Ok(())
}

fn run_analyze(args: &[String]) -> Result<(), String> {
    let mut archive: Option<PathBuf> = None;
    let mut panel: Option<PathBuf> = None;
    let mut json = false;
    let mut bounded = false;
    let mut threads = std::thread::available_parallelism().map_or(4, |n| n.get() as i32);

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--archive" | "--input" => archive = Some(PathBuf::from(next_val(args, &mut i, "--archive")?)),
            "--panel" | "--ref" => panel = Some(PathBuf::from(next_val(args, &mut i, "--panel")?)),
            "--json" => json = true,
            "--bounded" => bounded = true,
            "--threads" => { threads = next_val(args, &mut i, "--threads")?.parse::<i32>().map_err(|_| "Invalid --threads".to_string())?.max(1); }
            other => return Err(format!("Unknown analyze option: {other}")),
        }
        i += 1;
    }
    let arc = archive.ok_or("analyze requires --archive <path>")?;
    let rep = compute::analyze_archive(&arc, panel.as_deref(), threads, bounded)?;
    if json {
        println!("{}", rep.to_json());
    } else {
        print!("{}", rep.format_summary());
    }
    Ok(())
}

fn run_variants(args: &[String]) -> Result<(), String> {
    let mut archive: Option<PathBuf> = None;
    let mut panel: Option<PathBuf> = None;
    let mut out_vcf: Option<PathBuf> = None;
    let mut json = false;
    let mut bounded = false;
    let mut threads = std::thread::available_parallelism().map_or(4, |n| n.get() as i32);

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--archive" | "--input" => archive = Some(PathBuf::from(next_val(args, &mut i, "--archive")?)),
            "--panel" | "--ref" => panel = Some(PathBuf::from(next_val(args, &mut i, "--panel")?)),
            "--out" | "--output" => out_vcf = Some(PathBuf::from(next_val(args, &mut i, "--out")?)),
            "--json" => json = true,
            "--bounded" => bounded = true,
            "--threads" => { threads = next_val(args, &mut i, "--threads")?.parse::<i32>().map_err(|_| "Invalid --threads".to_string())?.max(1); }
            other => return Err(format!("Unknown variants option: {other}")),
        }
        i += 1;
    }
    let arc = archive.ok_or("variants requires --archive <path>")?;
    let res = compute::extract_variants_vcf(&arc, panel.as_deref(), out_vcf.as_deref(), threads, bounded)?;
    if json {
        println!(
            "{{\"total_variants\":{},\"snvs\":{},\"insertions\":{},\"deletions\":{},\"complex_indels\":{},\"literal_runs_skipped\":{},\"elapsed_seconds\":{:.6}}}",
            res.total_variants, res.snvs, res.insertions, res.deletions, res.complex_indels, res.literal_runs_skipped, res.seconds
        );
    } else if out_vcf.is_some() {
        println!(
            "VCF emitted: {} variants (SNV={}, INS={}, DEL={}, COMPLEX={}) in {:.4} s",
            res.total_variants, res.snvs, res.insertions, res.deletions, res.complex_indels, res.seconds
        );
    } else {
        print!("{}", res.vcf_text);
    }
    Ok(())
}

fn run_kmer(args: &[String]) -> Result<(), String> {
    let mut archive: Option<PathBuf> = None;
    let mut panel: Option<PathBuf> = None;
    let mut query: Option<String> = None;
    let mut both_strands = false;
    let mut json = false;
    let mut bounded = false;
    let mut threads = std::thread::available_parallelism().map_or(4, |n| n.get() as i32);

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--archive" | "--input" => archive = Some(PathBuf::from(next_val(args, &mut i, "--archive")?)),
            "--panel" | "--ref" => panel = Some(PathBuf::from(next_val(args, &mut i, "--panel")?)),
            "--query" | "-q" => query = Some(next_val(args, &mut i, "--query")?),
            "--both-strands" => both_strands = true,
            "--json" => json = true,
            "--bounded" => bounded = true,
            "--threads" => { threads = next_val(args, &mut i, "--threads")?.parse::<i32>().map_err(|_| "Invalid --threads".to_string())?.max(1); }
            other => return Err(format!("Unknown kmer option: {other}")),
        }
        i += 1;
    }
    let arc = archive.ok_or("kmer requires --archive <path> and --query <ACGT>")?;
    let q = query.ok_or("kmer requires --query <ACGT>")?;
    let rep = compute::search_kmers(&arc, panel.as_deref(), &q, both_strands, threads, bounded)?;
    if json {
        println!(
            "{{\"query_len\":{},\"hit_count\":{},\"exact_zero_target_spans\":{},\"patched_or_local_spans\":{},\"elapsed_seconds\":{:.6}}}",
            rep.query_len, rep.hits.len(), rep.exact_zero_target_spans, rep.patched_or_local_spans, rep.seconds
        );
    } else {
        println!(
            "Found {} exact hits for {}-mer (zero-target-spans={}, patched-spans={}, {:.4} s):",
            rep.hits.len(), rep.query_len, rep.exact_zero_target_spans, rep.patched_or_local_spans, rep.seconds
        );
        for h in rep.hits.iter().take(50) {
            println!("  unit={}\tpos={}\tstrand={}", h.unit_index, h.tgt_pos, if h.strand == 0 { '+' } else { '-' });
        }
    }
    Ok(())
}

fn run_query(args: &[String]) -> Result<(), String> {
    let mut archive: Option<PathBuf> = None;
    let mut panel: Option<PathBuf> = None;
    let mut unit = 0usize;
    let mut start: Option<u64> = None;
    let mut end: Option<u64> = None;
    let mut bounded = false;
    let mut threads = std::thread::available_parallelism().map_or(4, |n| n.get() as i32);

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--archive" | "--input" => archive = Some(PathBuf::from(next_val(args, &mut i, "--archive")?)),
            "--panel" | "--ref" => panel = Some(PathBuf::from(next_val(args, &mut i, "--panel")?)),
            "--unit" => { unit = next_val(args, &mut i, "--unit")?.parse().map_err(|_| "Invalid --unit".to_string())?; }
            "--start" => { start = Some(next_val(args, &mut i, "--start")?.parse().map_err(|_| "Invalid --start".to_string())?); }
            "--end" => { end = Some(next_val(args, &mut i, "--end")?.parse().map_err(|_| "Invalid --end".to_string())?); }
            "--bounded" => bounded = true,
            "--threads" => { threads = next_val(args, &mut i, "--threads")?.parse::<i32>().map_err(|_| "Invalid --threads".to_string())?.max(1); }
            other => return Err(format!("Unknown query option: {other}")),
        }
        i += 1;
    }
    let arc = archive.ok_or("query requires --archive <path>")?;
    let s = start.ok_or("query requires --start <pos>")?;
    let e = end.ok_or("query requires --end <pos>")?;
    let sub = compute::query_slice(&arc, panel.as_deref(), unit, s, e, threads, bounded)?;
    println!("{sub}");
    Ok(())
}

fn run_inspect(args: &[String]) -> Result<(), String> {
    let mut archive: Option<PathBuf> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--archive" | "--input" => archive = Some(PathBuf::from(next_val(args, &mut i, "--archive")?)),
            other => return Err(format!("Unknown inspect option: {other}")),
        }
        i += 1;
    }
    let arc = archive.ok_or("inspect requires --archive <path>")?;
    let info = hgs2::inspect_file(&arc)?;
    println!(
        "{{\"mode\":{},\"raw_flat_acgt\":{},\"raw_rna_u\":{},\"raw_bytes\":{},\"crc32\":{},\"total_records\":{},\"num_units\":{},\"total_novelty_bases\":{}}}",
        info.mode, info.raw_flat_acgt, info.raw_rna_u, info.raw_bytes, info.crc32, info.total_records, info.num_units, info.total_novelty_bases
    );
    Ok(())
}
