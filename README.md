# reaper
The REAlign error PronE Reads variant caller.

Reaper is a prototype SNV caller for long error prone reads such as Pacific Biosciences SMRT.

## dependencies
* rust 1.26.2 (see installation)
* various rust dependencies (automatically managed by cargo)

## installation
First, install the Rust programming language. Reaper requires Rust 1.26.2 which comes with the cargo package manager. You can install rust with the following command:
```
curl https://sh.rustup.rs -sSf | sh
```
To build reaper, type:
```
cargo build --release
```
This will compile the code with optimizations in release mode, and the binary will be
in ```target/release/reaper```. It will automatically install Rust dependencies for reaper,
such rust-bio and rust-htslib.

To run unit tests, type:
```
cargo test
```

Usage:
```
$ ./target/release/reaper --help
```
