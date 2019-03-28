// build.rs

extern crate cc;
extern crate cmake;

use cmake::Config;
use std::env;

fn main() {
    cc::Build::new()
        .flag_if_supported("-O3")
        .flag_if_supported("-Wall")
        .file("src/hapcut2/common.c")
        .file("src/hapcut2/fragmatrix.c")
        .file("src/hapcut2/readinputbuffers.c")
        .file("src/hapcut2/pointerheap.c")
        .file("src/hapcut2/hapcut2.c")
        .compile("hapcut2");

    let dst = Config::new("src/poa/spoa")
           .define("CMAKE_BUILD_TYPE","Release")
           .build();

    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=static=spoa");

    let out_dir = env::var("OUT_DIR").unwrap();
    println!("cargo:rustc-flags=-L {}/lib64/", &out_dir);

    cc::Build::new()
        .cpp(true)
        .shared_flag(false)
        .static_flag(true)
        .flag_if_supported("-O3")
        .flag_if_supported("-D_GNU_SOURCE")
        .flag_if_supported("-Wall")
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-Isrc/poa/spoa/include")
        .flag_if_supported(&format!("-L{}/lib64", &out_dir))
        .flag_if_supported("-lspoa")
        .file("src/poa/poa_func.cpp")
        .compile("poa_func");

    //println!("cargo:rustc-flags=-L /home/pedge/anaconda3/envs/tscc/lib");
}
