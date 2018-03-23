# reaper
The REAlign error PronE Reads variant caller.

Reaper is a prototype SNV caller for long error prone reads such as Pacific Biosciences SMRT.

## dependencies
* gcc/g++ >= 4.8
* cmake
* libclang
* dependencies of HTSlib (TODO)
* spoa (see installation)
* rust 1.4 (see installation)
* various rust dependencies (automatically managed by cargo)

## pre-installation

TODO describe the various dependencies for HTSlib and provide example command to install them with apt-get on ubuntu.

## installation
First, clone the repository with --recursive (this downloads Reaper and the SPOA library dependency):
```
git clone --recursive https://github.com/pjedge/reaper
cd reaper
```
Make sure your C and C++ compilers are properly specified:
```
export CC=gcc
export CXX=g++
```
Next, build the SPOA library (requires cmake) with position independent code on: 
```
mkdir src/poa/spoa/build
cd src/poa/spoa/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=TRUE ..
make
cd ../../../../
```
Next, install the Rust programming language. Reaper requires Rust 1.4 which comes with the cargo package manager. You can install rust with the following command:
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

## troubleshooting
### linker errors
For example:
```
error: linking with `cc` failed: exit code: 1
...
...
...
= note: src/poa/spoa/build/lib//libspoa.a(graph.cpp.o): In function `spoa::Node::coverage() const':
          graph.cpp:(.text+0x1ce5): undefined reference to `std::__detail::_Prime_rehash_policy::_M_next_bkt(unsigned long) const'
...
...
...
```
Your system may have multiple versions of your linker that are causing a conflict. Rustc may be calling to a different or old version of the linker. In this case, specify the linker (in linux, gcc) as follows:
```
rustc -vV
```
Note the build target after "host: ", i.e. "x86_64-unknown-linux-gnu".
```
mkdir .cargo
nano .cargo/config
```
edit the config file to have these contents:
```
[target.<target-name>]
linker = "</path/to/linker>"
```
for example,
```
[target.x86_64-unknown-linux-gnu]
linker = "/opt/gnu/gcc/bin/gcc"
```
then,
```
cargo clean
cargo build --release
```
