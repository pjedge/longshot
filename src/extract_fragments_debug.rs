use bio::io::fasta;
use extract_fragments::*;
use realignment::*;
use rust_htslib::bam;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use util::*;
use variants_and_fragments::*;

static IGNORE_INDEL_ONLY_CLUSTERS: bool = false;

//************************************************************************************************
// BEGINNING OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************

// This function (create_augmented_cigarlist) is a modification of Rust-Htslib function rust_htslib::bam::record::CigarStringView::read_pos
// https://github.com/rust-bio/rust-htslib/blob/master/src/bam/record.rs
// Rust-Htslib license is copied here as per its terms:
//The MIT License (MIT)
//
//Copyright (c) 2016 Johannes Köster, the Rust-Htslib team.

//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

pub fn get_fragment_anchors(
    bam_record: &Record,
    cigarpos_list: &Vec<CigarPos>,
    vars: Vec<Var>,
    ref_seq: &Vec<char>,
    target_names: &Vec<String>,
    extract_params: ExtractFragmentParameters,
    _align_params: AlignmentParameters,
    _old_frag: Option<Fragment>,
) -> Vec<(AnchorPositions, Vec<Var>)> {
    // TODO assert that every single variant in vars is on the same chromosome
    //let id: String = u8_to_string(bam_record.qname());

    if bam_record.is_quality_check_failed()
        || bam_record.is_duplicate()
        || bam_record.is_secondary()
        || bam_record.is_unmapped()
        || bam_record.mapq() < extract_params.min_mapq
        || bam_record.is_supplementary()
    {
        return vec![];
    }

    let read_seq: Vec<char> = dna_vec(&bam_record.seq().as_bytes());
    let mut cluster_lst: Vec<(AnchorPositions, Vec<Var>)> = vec![];
    let mut var_anchor_lst: Vec<(Var, AnchorPositions)> = vec![];

    // populate a list with tuples of each variant, and anchor sequences for its alignment
    for ref var in vars {
        let var_interval = GenomicInterval {
            tid: var.tid as u32,
            chrom: var.chrom.clone(),
            start_pos: var.pos0 as u32,
            end_pos: var.pos0 as u32,
        };
        match find_anchors(
            &bam_record,
            &cigarpos_list,
            var_interval,
            &ref_seq,
            &read_seq,
            &target_names,
            extract_params,
        )
        .expect("CIGAR or Anchor Error while finding anchor sequences.")
        {
            Some(anchors) => {
                var_anchor_lst.push((var.clone(), anchors));
            }
            _ => {}
        };
    }

    // now that we have anchors for each var the read covers,
    // group the variants into clusters to align together if adjacent anchors overlap
    // we generate anchors for the whole cluster by taking the first-left and last-right anchor pos of the cluster
    {
        // populate cluster_lst with variant clusters
        let mut var_cluster: Vec<Var> = vec![];
        let mut var_anchors: Vec<AnchorPositions> = vec![];
        let mut l;

        // generate clusters of SNVs that should be considered
        for (var, anc) in var_anchor_lst {
            l = var_cluster.len();
            if l == 0
                || (anc.left_anchor_ref < var_anchors[l - 1].right_anchor_ref
                    && l < extract_params.short_hap_max_snvs)
            {
                var_cluster.push(var);
                var_anchors.push(anc);
            } else {
                // sequence anchor that covers the whole cluster of variants
                let combined_anchor = AnchorPositions {
                    left_anchor_ref: var_anchors[0].left_anchor_ref,
                    right_anchor_ref: var_anchors[l - 1].right_anchor_ref,
                    left_anchor_read: var_anchors[0].left_anchor_read,
                    right_anchor_read: var_anchors[l - 1].right_anchor_read,
                };

                cluster_lst.push((combined_anchor, var_cluster.clone()));

                var_cluster.clear();
                var_anchors.clear();
                var_cluster.push(var);
                var_anchors.push(anc);
            }
        }

        l = var_cluster.len();

        if l > 0 {
            // sequence anchor that covers the whole cluster of variants
            let combined_anchor = AnchorPositions {
                left_anchor_ref: var_anchors[0].left_anchor_ref,
                right_anchor_ref: var_anchors[l - 1].right_anchor_ref,
                left_anchor_read: var_anchors[0].left_anchor_read,
                right_anchor_read: var_anchors[l - 1].right_anchor_read,
            };

            cluster_lst.push((combined_anchor, var_cluster.clone()));
        }
    }

    cluster_lst
}

pub fn extract_fragments_debug(
    bam_file: &String,
    fastafile_name: &String,
    varlist: &VarList,
    interval: &Option<GenomicInterval>,
    extract_params: ExtractFragmentParameters,
    align_params: AlignmentParameters,
    old_flist: Option<Vec<Fragment>>,
) -> () {
    let t_names = parse_target_names(&bam_file);

    let mut prev_tid = 4294967295; // huge value so that tid != prev_tid on first iter
    let mut fasta = fasta::IndexedReader::from_file(fastafile_name).unwrap();
    let mut ref_seq: Vec<char> = vec![];

    let mut total_cluster_lst = vec![];

    // TODO: this uses a lot of duplicate code, need to figure out a better solution.
    let mut complete = 0;

    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_file, interval);
    let mut bam_ix = bam::IndexedReader::from_path(bam_file).unwrap();

    for iv in interval_lst {
        bam_ix
            .fetch(iv.tid, iv.start_pos, iv.end_pos + 1)
            .ok()
            .expect("Error seeking BAM file while extracting fragments.");

        for (i, r) in bam_ix.records().enumerate() {
            let record = r.unwrap();

            if record.is_quality_check_failed()
                || record.is_duplicate()
                || record.is_secondary()
                || record.is_unmapped()
                || record.mapq() < extract_params.min_mapq
                || record.is_supplementary()
            {
                continue;
            }

            let old_frag: Option<Fragment> = match &old_flist {
                &Some(ref fl) => {
                    assert_eq!(fl[i].id, u8_to_string(record.qname()));
                    Some(fl[i].clone())
                }
                &None => None,
            };

            let tid: usize = record.tid() as usize;
            let chrom: String = t_names[record.tid() as usize].clone();

            if tid != prev_tid {
                let mut ref_seq_u8: Vec<u8> = vec![];
                fasta
                    .read_all(&chrom, &mut ref_seq_u8)
                    .expect("Failed to read fasta sequence record.");
                ref_seq = dna_vec(&ref_seq_u8);
            }

            let start_pos = record.pos();
            let end_pos = record
                .cigar()
                .end_pos()
                .expect("Error while accessing CIGAR end position")
                - 1;

            let bam_cig: CigarStringView = record.cigar();
            let cigarpos_list: Vec<CigarPos> =
                create_augmented_cigarlist(start_pos as u32, &bam_cig)
                    .expect("Error creating augmented cigarlist.");

            let interval = GenomicInterval {
                tid: tid as u32,
                chrom: chrom,
                start_pos: start_pos as u32,
                end_pos: end_pos as u32,
            };

            // get the list of variants that overlap this read
            let read_vars = varlist.get_variants_range(interval);

            // print the percentage of variants processed every 10%
            if read_vars.len() > 0
                && ((read_vars[0].ix as f64 / varlist.lst.len() as f64) * 10.0) as usize > complete
            {
                complete = ((read_vars[0].ix as f64 / varlist.lst.len() as f64) * 10.0) as usize;
                if complete < 10 {
                    eprintln!(
                        "{}    {}% of variants processed...",
                        print_time(),
                        complete * 10
                    );
                }
            }

            let mut cluster_lst = get_fragment_anchors(
                &record,
                &cigarpos_list,
                read_vars,
                &ref_seq,
                &t_names,
                extract_params,
                align_params,
                old_frag,
            );

            total_cluster_lst.append(&mut cluster_lst);

            prev_tid = tid;
        }
    }
    eprintln!("{}    100% of variants processed.", print_time());

    let mut total_ops = 0;
    let mut total_ref_len = 0;
    let mut total_read_len = 0;
    let mut total_variants = 0;
    let region_width = 500000;

    let max_pos = total_cluster_lst.iter().map(|x| x.1[0].pos0).max().unwrap();
    let max_ix = max_pos / region_width + 1;

    let mut region_op_counts = vec![0; max_ix];
    let mut region_ref_len_sum = vec![0; max_ix];
    let mut region_read_len_sum = vec![0; max_ix];
    let mut region_variants = vec![0; max_ix];
    let mut region_cluster_count = vec![0; max_ix];

    for &(ref anchor, ref var_cluster) in &total_cluster_lst {
        // estimate the number of operations used in alignment
        let ref_len = anchor.right_anchor_ref - anchor.left_anchor_ref + 1;
        let read_len = anchor.right_anchor_read - anchor.left_anchor_read + 1;
        let len_diff = ((read_len as i32) - (ref_len as i32)).abs() as usize;
        let band_width = extract_params.band_width + len_diff;
        let n_ops = read_len as usize * band_width;

        // multiply by the number of haplotypes being aligned to
        let n_haps = (2 as usize).pow(var_cluster.len() as u32);

        // update the total counts
        total_ops += n_haps * n_ops;
        total_ref_len += ref_len;
        total_read_len += read_len;
        total_variants += var_cluster.len();

        // get the index corresponding to this var cluster's "region"
        let region_ix = var_cluster[0].pos0 / region_width;

        // update the counts for specific regions
        region_op_counts[region_ix] += n_haps * n_ops;
        region_ref_len_sum[region_ix] += ref_len;
        region_read_len_sum[region_ix] += read_len;
        region_variants[region_ix] += var_cluster.len();
        region_cluster_count[region_ix] += 1;
    }
    eprintln!("\n----------------------------------------------------------------");
    eprintln!("RESULTS FOR ALLELE REALIGNMENT DEBUG ANALYSIS:");

    let total_mean_read_window_len = total_read_len as f64 / (total_cluster_lst.len() as f64);
    let total_mean_ref_window_len = total_ref_len as f64 / (total_cluster_lst.len() as f64);
    let mean_variants_per_cluster = total_variants as f64 / (total_cluster_lst.len() as f64);

    eprintln!("Total realignment operations:  {}", total_ops);
    eprintln!(
        "Total mean read window length: {:.2}",
        total_mean_read_window_len
    );
    eprintln!(
        "Total mean ref window length:  {:.2}",
        total_mean_ref_window_len
    );
    eprintln!(
        "Mean variants per cluster:     {:.2}",
        mean_variants_per_cluster
    );

    eprintln!("\n----------------------------------------------------------------");
    eprintln!("RESULTS PER REGION:");
    eprintln!("REGION\tMEAN_VARIANTS_PER_CLUSTER\tMEAN_READ_WINDOW_LEN\tMEAN_REF_WINDOW_LEN\tOPS_COUNT\tOPS_ENRICHMENT");

    let mean_ops_per_region = total_ops as f64 / (region_op_counts.len() as f64);

    for i in 0..max_ix {
        let chrom = interval.clone().unwrap().chrom;
        let region_start = (i * region_width) as usize;
        let region_end = ((i + 1) * region_width) - 1 as usize;
        let ops_count = region_op_counts[i];
        let ops_enrichment = region_op_counts[i] as f64 / mean_ops_per_region;
        let mean_variants_per_cluster =
            region_variants[i] as f64 / (region_cluster_count[i] as f64);
        let mean_read_window_len = region_read_len_sum[i] as f64 / (region_cluster_count[i] as f64);
        let mean_ref_window_len = region_ref_len_sum[i] as f64 / (region_cluster_count[i] as f64);
        eprintln!(
            "{}:{}-{}\t{:.2}\t{:.2}\t{:.2}\t{}\t{:.2}",
            chrom,
            region_start,
            region_end,
            mean_variants_per_cluster,
            mean_read_window_len,
            mean_ref_window_len,
            ops_count,
            ops_enrichment
        );
    }
}

//************************************************************************************************
// END OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************
