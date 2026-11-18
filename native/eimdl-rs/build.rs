use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

fn tool_exists(cmd: &str, arg: &str) -> bool {
    Command::new(cmd)
        .arg(arg)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

fn run_checked(cmd: &mut Command, label: &str) {
    let status = cmd
        .status()
        .unwrap_or_else(|e| panic!("Failed to launch {label}: {e}"));
    if !status.success() {
        panic!("{label} failed with exit status {status}");
    }
}

fn pick_std_flag(
    tool: &str,
    prefix: &[&str],
    is_cpp: bool,
    out_dir: &Path,
    candidates: &[&str],
) -> &'static str {
    let ext = if is_cpp { "cpp" } else { "c" };
    let probe_src = out_dir.join(format!("_std_probe.{ext}"));
    let probe_obj = out_dir.join(format!("_std_probe_{ext}.o"));
    let body = if is_cpp {
        "#include <bit>\nint _hgs_std_probe() { return std::countl_zero(1U); }\n"
    } else {
        "int _hgs_std_probe(void) { return 0; }\n"
    };
    let _ = std::fs::write(&probe_src, body);
    let mut chosen: Option<&'static str> = None;
    for &flag in candidates {
        let mut cmd = Command::new(tool);
        cmd.args(prefix);
        if is_cpp {
            cmd.args(["-x", "c++"]);
        } else {
            cmd.args(["-x", "c"]);
        }
        cmd.arg(flag).args(["-c"]).arg(&probe_src).arg("-o").arg(&probe_obj);
        if let Ok(out) = cmd.output() {
            if out.status.success() {
                chosen = match flag {
                    "-std=c23" => Some("-std=c23"),
                    "-std=c2x" => Some("-std=c2x"),
                    "-std=c++23" => Some("-std=c++23"),
                    "-std=c++2b" => Some("-std=c++2b"),
                    _ => None,
                };
                if chosen.is_some() {
                    break;
                }
            }
        }
    }
    let _ = std::fs::remove_file(&probe_src);
    let _ = std::fs::remove_file(&probe_obj);
    chosen.unwrap_or_else(|| {
        panic!(
            "Compiler '{tool}' supports none of the required modern standards: {}",
            candidates.join(", ")
        )
    })
}

fn collect_zig_cxx_libs(zig_exe: &str, arch_target: &str, cxx_std: &str, out_dir: &Path) -> Vec<(&'static str, PathBuf)> {
    // Ask `zig c++ -v -shared` to resolve the exact C++ runtime archives for `arch_target`
    // rather than scanning the global zig cache (which may hold artifacts for other targets such as wasm32).
    let probe_cpp = out_dir.join("_zig_probe.cpp");
    let probe_dll = out_dir.join(if arch_target.contains("windows") { "_zig_probe.dll" } else { "_zig_probe.so" });
    let _ = std::fs::write(
        &probe_cpp,
        "#include <vector>\n#include <thread>\nextern \"C\" int _hgs_zprobe(){std::vector<int> v(1);return (int)v.size();}\n",
    );
    let mut found: Vec<(&'static str, PathBuf)> = Vec::new();
    if let Ok(out) = Command::new(zig_exe)
        .args(["c++", "-target", arch_target, cxx_std, "-fno-sanitize=undefined", "-v", "-shared", "-o"])
        .arg(&probe_dll)
        .arg(&probe_cpp)
        .output()
    {
        let combined = format!(
            "{}\n{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        let targets: [(&str, &[&str]); 3] = [
            ("c++", &["libc++.a", "c++.lib"]),
            ("c++abi", &["libc++abi.a", "c++abi.lib"]),
            ("unwind", &["libunwind.a", "unwind.lib"]),
        ];
        for (rust_name, candidates) in targets {
            'tok: for line in combined.lines() {
                for tok in line.split_whitespace() {
                    let clean = tok.trim_matches(|c| c == '"' || c == '\'' || c == ',');
                    let p = Path::new(clean);
                    if let Some(fname) = p.file_name().and_then(|s| s.to_str()) {
                        if candidates.iter().any(|c| fname.eq_ignore_ascii_case(c)) && p.exists() {
                            found.push((rust_name, p.to_path_buf()));
                            break 'tok;
                        }
                    }
                }
            }
        }
    }
    let _ = std::fs::remove_file(&probe_cpp);
    let _ = std::fs::remove_file(&probe_dll);
    if arch_target.contains("windows") {
        let _ = std::fs::remove_file(out_dir.join("_zig_probe.pdb"));
        let _ = std::fs::remove_file(out_dir.join("_zig_probe.lib"));
    }
    found
}

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let native_dir = manifest_dir.parent().unwrap().to_path_buf();
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();

    let c_sources = ["hgs_simd.c", "hgs_deflate.c", "hgs_ans.c"];
    let c_headers = ["hgs_simd.h", "hgs_deflate.h", "hgs_ans.h"];
    let cpp_src = native_dir.join("cpp_engine").join("hgs_cpp_engine.cpp");
    let c_inc = native_dir.join("c_kernels");
    let cpp_inc = native_dir.join("cpp_engine");

    for src in c_sources {
        println!("cargo:rerun-if-changed={}", c_inc.join(src).display());
    }
    for hdr in c_headers {
        println!("cargo:rerun-if-changed={}", c_inc.join(hdr).display());
    }
    println!("cargo:rerun-if-changed={}", cpp_src.display());
    println!("cargo:rerun-if-changed={}", cpp_inc.join("hgs_cpp_engine.hpp").display());

    let opt_flag = if env::var("PROFILE").as_deref() == Ok("debug") { "-O2" } else { "-O3" };
    let c_objs: Vec<PathBuf> = c_sources
        .iter()
        .map(|s| out_dir.join(s.replace(".c", ".o")))
        .collect();
    let cpp_obj = out_dir.join("hgs_cpp_engine.o");
    let lib_path = out_dir.join("libhgs_native.a");

    let env_cc = env::var("CC")
        .ok()
        .filter(|s| !s.trim().is_empty() && tool_exists(s, "--version"));
    let env_cxx = env::var("CXX")
        .ok()
        .filter(|s| !s.trim().is_empty() && tool_exists(s, "--version"));
    let has_env_cc = env_cc.is_some() && env_cxx.is_some();
    let has_gcc = tool_exists("gcc", "--version") && tool_exists("g++", "--version");
    let has_clang = tool_exists("clang", "--version") && tool_exists("clang++", "--version");
    let has_zig = tool_exists("zig", "version");

    if has_env_cc || has_gcc || (target_os != "windows" && has_clang) {
        let prefer_clang = target_os == "macos" && has_clang;
        let cc = env_cc.unwrap_or_else(|| {
            if prefer_clang || !has_gcc { "clang".into() } else { "gcc".into() }
        });
        let cxx = env_cxx.unwrap_or_else(|| {
            if prefer_clang || !has_gcc { "clang++".into() } else { "g++".into() }
        });
        let ar = env::var("AR")
            .ok()
            .filter(|s| !s.trim().is_empty())
            .unwrap_or_else(|| "ar".into());
        let c_std = pick_std_flag(&cc, &[], false, &out_dir, &["-std=c23", "-std=c2x"]);
        let cxx_std = pick_std_flag(&cxx, &[], true, &out_dir, &["-std=c++23", "-std=c++2b"]);

        for (src, obj) in c_sources.iter().zip(c_objs.iter()) {
            run_checked(
                Command::new(&cc)
                    .args([c_std, opt_flag, "-Wall", "-Wextra", "-fPIC", "-I"])
                    .arg(&c_inc)
                    .args(["-c"])
                    .arg(c_inc.join(src))
                    .arg("-o")
                    .arg(obj),
                &format!("compile {src}"),
            );
        }
        run_checked(
            Command::new(&cxx)
                .args([cxx_std, opt_flag, "-Wall", "-Wextra", "-fPIC", "-I"])
                .arg(&c_inc)
                .arg("-I")
                .arg(&cpp_inc)
                .args(["-c"])
                .arg(&cpp_src)
                .arg("-o")
                .arg(&cpp_obj),
            "compile hgs_cpp_engine.cpp",
        );
        let _ = std::fs::remove_file(&lib_path);
        let mut ar_cmd = Command::new(&ar);
        ar_cmd.args(["rcs"]).arg(&lib_path);
        for obj in &c_objs {
            ar_cmd.arg(obj);
        }
        ar_cmd.arg(&cpp_obj);
        run_checked(&mut ar_cmd, "archive libhgs_native.a");

        println!("cargo:rustc-link-search=native={}", out_dir.display());
        println!("cargo:rustc-link-lib=static=hgs_native");
        if target_os == "macos" || cxx.contains("clang") {
            println!("cargo:rustc-link-lib=c++");
        } else {
            println!("cargo:rustc-link-lib=stdc++");
        }
        if target_os == "linux" {
            println!("cargo:rustc-link-lib=pthread");
            println!("cargo:rustc-link-lib=m");
        }
    } else if has_zig {
        let zig_target = if target_os == "windows" && target_env == "gnu" {
            "x86_64-windows-gnu"
        } else if target_os == "linux" {
            "x86_64-linux-gnu"
        } else {
            "native"
        };
        let c_std = pick_std_flag("zig", &["cc", "-target", zig_target], false, &out_dir, &["-std=c23", "-std=c2x"]);
        let cxx_std = pick_std_flag("zig", &["c++", "-target", zig_target], true, &out_dir, &["-std=c++23", "-std=c++2b"]);
        if target_os == "windows" {
            for (src, obj) in c_sources.iter().zip(c_objs.iter()) {
                run_checked(
                    Command::new("zig")
                        .args(["cc", "-target", zig_target, c_std, opt_flag, "-fno-sanitize=undefined", "-Wall", "-Wextra", "-I"])
                        .arg(&c_inc)
                        .args(["-c"])
                        .arg(c_inc.join(src))
                        .arg("-o")
                        .arg(obj),
                    &format!("zig cc {src}"),
                );
            }
            run_checked(
                Command::new("zig")
                    .args(["c++", "-target", zig_target, cxx_std, opt_flag, "-fno-sanitize=undefined", "-Wall", "-Wextra", "-I"])
                    .arg(&c_inc)
                    .arg("-I")
                    .arg(&cpp_inc)
                    .args(["-c"])
                    .arg(&cpp_src)
                    .arg("-o")
                    .arg(&cpp_obj),
                "zig c++ hgs_cpp_engine.cpp",
            );
            let dll_path = out_dir.join("hgs_native.dll");
            let imp_path = out_dir.join("libhgs_native.dll.a");
            let mut imp_arg = std::ffi::OsString::from("-Wl,--out-implib,");
            imp_arg.push(&imp_path);
            let mut link_cmd = Command::new("zig");
            link_cmd.args(["c++", "-target", zig_target, cxx_std, opt_flag, "-fno-sanitize=undefined", "-shared"]);
            for obj in &c_objs {
                link_cmd.arg(obj);
            }
            link_cmd.arg(&cpp_obj).arg("-Wl,--export-all-symbols").arg(&imp_arg).arg("-o").arg(&dll_path);
            run_checked(&mut link_cmd, "zig c++ link hgs_native.dll");
            if let Some(build_dir) = out_dir.ancestors().find(|p| p.file_name().and_then(|s| s.to_str()) == Some("build")) {
                if let Some(profile_dir) = build_dir.parent() {
                    let _ = std::fs::copy(&dll_path, profile_dir.join("hgs_native.dll"));
                    let deps_dir = profile_dir.join("deps");
                    let _ = std::fs::create_dir_all(&deps_dir);
                    let _ = std::fs::copy(&dll_path, deps_dir.join("hgs_native.dll"));
                }
            }
            println!("cargo:rustc-link-search=native={}", out_dir.display());
            println!("cargo:rustc-link-lib=dylib=hgs_native");
        } else {
            for (src, obj) in c_sources.iter().zip(c_objs.iter()) {
                run_checked(
                    Command::new("zig")
                        .args(["cc", "-target", zig_target, c_std, opt_flag, "-fno-sanitize=undefined", "-fPIC", "-Wall", "-Wextra", "-I"])
                        .arg(&c_inc)
                        .args(["-c"])
                        .arg(c_inc.join(src))
                        .arg("-o")
                        .arg(obj),
                    &format!("zig cc {src}"),
                );
            }
            run_checked(
                Command::new("zig")
                    .args(["c++", "-target", zig_target, cxx_std, opt_flag, "-fno-sanitize=undefined", "-fPIC", "-Wall", "-Wextra", "-I"])
                    .arg(&c_inc)
                    .arg("-I")
                    .arg(&cpp_inc)
                    .args(["-c"])
                    .arg(&cpp_src)
                    .arg("-o")
                    .arg(&cpp_obj),
                "zig c++ hgs_cpp_engine.cpp",
            );
            let _ = std::fs::remove_file(&lib_path);
            let mut ar_cmd = Command::new("zig");
            ar_cmd.args(["ar", "rcs"]).arg(&lib_path);
            for obj in &c_objs {
                ar_cmd.arg(obj);
            }
            ar_cmd.arg(&cpp_obj);
            run_checked(&mut ar_cmd, "zig ar libhgs_native.a");

            let cxx_libs = collect_zig_cxx_libs("zig", zig_target, cxx_std, &out_dir);
            println!("cargo:rustc-link-search=native={}", out_dir.display());
            println!("cargo:rustc-link-lib=static:-bundle=hgs_native");
            for (rust_name, src_path) in cxx_libs {
                if rust_name == "unwind" { continue; }
                let dst = out_dir.join(format!("lib{rust_name}.a"));
                std::fs::copy(&src_path, &dst).unwrap_or_else(|e| {
                    panic!("Failed to stage {} -> {}: {e}", src_path.display(), dst.display())
                });
                println!("cargo:rustc-link-lib=static:-bundle={rust_name}");
            }
            if target_os == "linux" {
                println!("cargo:rustc-link-lib=pthread");
                println!("cargo:rustc-link-lib=m");
            }
        }
    } else {
        panic!("No C23/C++23 compiler found (tried gcc/g++, clang/clang++, and zig)");
    }
}
