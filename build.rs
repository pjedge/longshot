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
}
