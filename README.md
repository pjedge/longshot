# longshot

Longshot is a diploid SNV caller for long error prone reads such as Pacific Biosciences SMRT.

THIS README IS A WORK IN PROGRESS AND MORE DETAILED DOCUMENTATION WILL BE AVAILABLE SOON

## supported operating systems
Longshot has been tested using Ubuntu 16.04, CentOS 6.6, and Manjaro Linux 17.1.11.
It should work on any linux-based system that has Rust and Cargo installed. 

## dependencies

* rust 1.26.2 (see installation)
* various rust dependencies (automatically managed by cargo)

## installation
First, install the Rust programming language. Longshot requires Rust 1.26.2 or higher which comes with the cargo package manager. You can install rust with the following command:
```
curl https://sh.rustup.rs -sSf | sh
```
To build Longshot, type:
```
cargo build --release
```
This will compile the code with optimizations in release mode, and the binary will be
in ```target/release/longshot```. It will automatically install Rust dependencies for reaper,
such rust-bio and rust-htslib. The build process should take around 3 minutes on a typical desktop computer.

To run unit tests, type:
```
cargo test
```

Usage:
```
$ ./target/release/longshot --help
```

## execution on an example dataset
The directory ```example_data``` contains a simulated toy dataset that can be used to test out Longshot:
- Reference genome containing 3 contigs each with length 200 kb (```example_data/genome.fa```)
- 30x coverage simulated pacbio reads generated using [SimLoRD](https://bitbucket.org/genomeinformatics/simlord/) (```example_data/pacbio_reads_30x.bam```)
- The 714 "true" variants for validation (```example_data/ground_truth_variants.vcf```)

Run Longshot on the example data as so:
```
./target/release/longshot --bam example_data/pacbio_reads_30x.bam --ref example_data/genome.fa --out example_data/longshot_output.vcf
```

Execution should take around 30 seconds on a typical desktop machine. The output can be compared to ```ground_truth_variants.vcf``` for accuracy.

## installation troubleshooting
### linker errors
For example:
```
error: linking with `cc` failed: exit code: 1
...
...
...
= note: Non-UTF-8 output: /usr/bin/ld: /home/pedge/temp/longshot/target/release/build/longshot-347f3774e75b380c/out/libhapcut2.a(common.o)(.text.fprintf_time+0x81): unresolvable H\x89\\$\xe8H\x89l$\xf0H\x89\xf3L\x89d$\xf8H\x83\xec\x18H\x8bG\x10H\x89\xfdI\x89\xd4H\x89\xd6H\x8b;\xffPxH\x8bE\x10I\x8dt$\x08H\x8b{\x08\xffPxH\x8bE\x10H\x8b{\x10I\x8dt$\x10H\x8b\x1c$H\x8bl$\x08L\x8bd$\x10H\x8b@xH\x83\xc4\x18\xff\xe0f\x90H\x89\\$\xe8H\x89l$\xf0H\x89\xfbL\x89d$\xf8H\x83\xec\x18H\x8bG\x10I\x89\xd4H\x89\xf5H\x89\xf7\xffPhI\x89\x04$H\x8bC\x10H\x8d}\x08\xffPhH\x8b\x1c$I\x89D$\x08H\x8bl$\x08L\x8bd$\x10H\x83\xc4\x18\xc3\x0f\x1f relocation against symbol `time@@GLIBC_2.2.5\'\n/usr/bin/ld: BFD version 2.20.51.0.2-5.42.el6 20100205 internal error, aborting at reloc.c line 443 in bfd_get_reloc_size\n\n/usr/bin/ld: Please report this bug.\n\ncollect2: ld returned 1 exit status\n
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
