// build.rs

extern crate cc;

fn main() {
    cc::Build::new()
        .flag_if_supported("-O3")
        .flag_if_supported("-D_GNU_SOURCE")
        .flag_if_supported("-Wall")
        .file("src/c/common.c")
        .file("src/c/fragmatrix.c")
        .file("src/c/readinputbuffers.c")
        .file("src/c/pointerheap.c")
        .file("src/c/hapcut2.c")
        .compile("hapcut2");

    cc::Build::new()
        .cpp(true)
        .static_flag(false)
        .flag_if_supported("-O3")
        .flag_if_supported("-D_GNU_SOURCE")
        .flag_if_supported("-Wall")
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-Isrc/poa/spoa/include/")
        .flag_if_supported("-Lsrc/poa/spoa/build/lib/")
        .flag_if_supported("-lspoa")
        .file("src/poa/poa_func.cpp")
        .compile("poa_func");

        println!("cargo:rustc-flags=-L src/poa/spoa/build/lib/ -l spoa");
        //println!("cargo:rustc-flags=-l dylib=stdc++ -lspoa");
    /*
    cc::Build::new()
        .flag_if_supported("-O3")
        .flag_if_supported("-D_GNU_SOURCE")
        .flag_if_supported("-Wall")
        .file("src/poa/poa_func.o")
        .file("src/poa/poa.c")
        .compile("poa");
    */

}
// .include("src/poa/spoa/include/")
// .flag_if_supported("-Lsrc/poa/spoa/build/lib/")
// .flag_if_supported("-lspoa")
