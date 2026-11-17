use std::env;
use std::fs;
use std::path::PathBuf;
use std::process;

use eimdl_native::ffi;
use eimdl_native::hgs2::{self, CompressOptions, MODE_COLLECTION, MODE_SINGLE};

fn print_usage() {
    eprintln!(
        "EIMDL / HyperGraSS v2.0.0 Native Engine (C23 AVX2/SWAR + C++23 64-bit DAG/CRM + Rust 2024 Streaming)\n\
         Usage:\n\
           eimdl --version\n\
           eimdl compress --target <file> --out <archive> [--panel <dir|file>] [--exclude <name>] [--panel-list <file>] [--panel-limit <N>] [--collection] [--single] [--bounded] [--threads <N>]\n\
           eimdl decompress --archive <archive> [--panel <dir|file>] [--bounded] [--threads <N>] --out <file>\n\
           eimdl inspect --archive <archive>"
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        print_usage();
        process::exit(1);
    }
    match args[1].as_str() {
        "--version" | "-V" | "version" => {
            println!("hypergrass-eimdl-native 2.0.0 (avx2={})", ffi::has_avx2());
        }
        "compress" => {
            if let Err(e) = run_compress(&args[2..]) {
                eprintln!("Error: {e}");
                process::exit(1);
            }
        }
        "decompress" => {
            if let Err(e) = run_decompress(&args[2..]) {
                eprintln!("Error: {e}");
                process::exit(1);
            }
        }
        "inspect" => {
            if let Err(e) = run_inspect(&args[2..]) {
                eprintln!("Error: {e}");
                process::exit(1);
            }
        }
        _ => {
            print_usage();
            process::exit(1);
        }
    }
}

fn next_val(args: &[String], i: &mut usize, flag: &str) -> Result<String, String> {
    *i += 1;
    if *i >= args.len() {
        return Err(format!("{flag} requires an argument"));
    }
    Ok(args[*i].clone())
}

fn run_compress(args: &[String]) -> Result<(), String> {
    let mut target: Option<PathBuf> = None;
    let mut out: Option<PathBuf> = None;
    let mut opts = CompressOptions::default();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--target" | "--input" => target = Some(PathBuf::from(next_val(args, &mut i, "--target")?)),
            "--out" | "--output" => out = Some(PathBuf::from(next_val(args, &mut i, "--out")?)),
            "--panel" | "--ref" => opts.reference_panel = Some(PathBuf::from(next_val(args, &mut i, "--panel")?)),
            "--exclude" => { opts.exclude_panel_names.insert(next_val(args, &mut i, "--exclude")?); }
            "--panel-list" => {
                let lp = PathBuf::from(next_val(args, &mut i, "--panel-list")?);
                let txt = fs::read_to_string(lp).map_err(|e| e.to_string())?;
                let list: Vec<String> = txt.lines().map(|l| l.trim().to_string()).filter(|l| !l.is_empty()).collect();
                opts.explicit_panel_names = Some(list);
            }
            "--panel-limit" => { opts.panel_limit = next_val(args, &mut i, "--panel-limit")?.parse().map_err(|_| "Invalid --panel-limit".to_string())?; }
            "--threads" => { opts.threads = next_val(args, &mut i, "--threads")?.parse::<i32>().map_err(|_| "Invalid --threads".to_string())?.max(1); }
            "--bounded" => opts.bounded = true,
            "--collection" => opts.mode = MODE_COLLECTION as i32,
            "--single" => { opts.single_genome = true; opts.mode = MODE_SINGLE as i32; }
            "--k" | "--step" => { let _ = next_val(args, &mut i, "--k")?; }
            "--no-approx" => opts.mosaic_flags.allow_approx = false,
            "--no-splice" => opts.mosaic_flags.allow_splice = false,
            "--no-v4" => opts.mosaic_flags.allow_kimura = false,
            "--no-wobble" => opts.mosaic_flags.allow_wobble = false,
            "--no-self" => opts.mosaic_flags.allow_self = false,
            "--no-tandem" => opts.mosaic_flags.allow_tandem = false,
            other => return Err(format!("Unknown compress option: {other}")),
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
    let secs = hgs2::decompress_file(&arc, &out_p, panel.as_deref(), threads, bounded)?;
    println!("{{\"decompress_seconds\":{secs:.6},\"engine\":\"native-c-cpp-rust\",\"verified_crc32\":true}}");
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
