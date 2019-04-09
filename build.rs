// build.rs

extern crate cc;

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
}
