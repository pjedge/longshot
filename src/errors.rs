
//! errors
//!
//! defines some custom errors to use with error-chain

error_chain! {
    errors {
        // BAM errors
        BamOpenError {
            description("Error opening BAM file.")
        }
        BamWriterOpenError(f: String) {
            description("Error opening BAM file for writing")
            display(x) -> ("{}: {}", x.description(), f)
        }
        BamRecordWriteError(qname: String) {
            description("Error writing record to bam file")
            display(x) -> ("{}: qname {}", x.description(), qname)
        }
        BamHeaderTargetLenAccessError {
            description("Error accessing target len for a contig in bam header.")
        }
        // Indexed BAM errors
        IndexedBamOpenError {
            description("Error opening indexed BAM file.")
        }
        IndexedBamReadError {
            description("Error reading indexed BAM file.")
        }
        IndexedBamFetchError {
            description("Error fetching region from indexed BAM file.")
        }
        IndexedBamRecordReadError {
            description("Error reading BAM record.")
        }
        IndexedBamPileupReadError {
            description("Error reading BAM pileup.")
        }
        IndexedBamPileupQueryPositionError {
            description("Error accessing query position for alignment in BAM pileup.")
        }
        // Indexed Fasta errors
        IndexedFastaOpenError {
            description("Error opening indexed FASTA file.")
        }
        IndexedFastaReadError {
            description("Error reading indexed FASTA file.")
        }
        // Indexed BAM errors
        BCFOpenError {
            description("Error opening BCF file.")
        }
        BCFReadError {
            description("Error reading BCF file.")
        }
        // CIGAR errors
        // derived from Rust-htslib errors defined with quick-error... https://github.com/rust-bio/rust-htslib/blob/master/src/bam/record.rs
        UnexpectedCigarOperation(msg: String) {
            description("Unexpected CIGAR operation")
            display(x) -> ("{}: {}", x.description(), msg)
        }
        UnsupportedCigarOperation(msg: String) {
            description("Unsupported CIGAR operation")
            display(x) -> ("{}: {}", x.description(), msg)
        }
        // Read realignment errors
        AnchorRangeOutsideRead {
            description("Attempted to find sequence anchors for a range completely outside of the sequence.")
        }
        // Gentotype priors errors
        InvalidTransitionBase(base: String) {
            description("Invalid base accessed from base transition hashmap")
            display(x) -> ("{}: {}", x.description(), base)
        }
        InvalidHaploidGenotype(refbase: char, altbase: char) {
            description("Invalid base accessed from base transition hashmap")
            display(x) -> ("{}: ref {}, alt {}", x.description(), refbase, altbase)
        }
        GenotypeNotInGenotypePriorsError(refbase: char, g0: char, g1: char) {
            description("Genotype not in genotype priors.")
            display(x) -> ("{}: ref {}, genotype ({},{})", x.description(), refbase, g0, g1)
        }
        // File IO errors
        FileWriteError (filename: String){
            description("Couldn't write to file")
            display(x) -> ("{}: {}", x.description(), filename)
        }
        CreateFileError (filename: String){
            description("Couldn't create file")
            display(x) -> ("{}: {}", x.description(), filename)
        }
        // File IO errors
        NoneError {
            description("Option was None.")
        }
    }
}
