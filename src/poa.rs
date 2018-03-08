// Copyright 2017 Brett Bowman
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Partial-Order Alignment for fast alignment and consensus of long
//! high-error reads (e.g. PacBio or Oxford Nanopore).  Represent the
//! multi-sequece alignment of many reads as a graph-structure, and use
//! Smith-Waterman alignments to the consensus-path to efficiently add
//! new sequences to the graph.  Because of this graph structure, both
//! traceback (i.e. consensus) is approximately O(N), and memory usage
//! is also O(N), where N is the number of nodes in the graph.
//!
//! For the original concept and theory, see:
//! * Lee, Christopher, Catherine Grasso, and Mark F. Sharlow. "Multiple sequence alignment using
//! partial order graphs." Bioinformatics 18.3 (2002): 452-464.
//! * Lee, Christopher. "Generating consensus sequences from partial order multiple sequence
//! alignment graphs." Bioinformatics 19.8 (2003): 999-1008.
//!
//! For a reference implementation that inspired this code, see poapy:
//! https://github.com/ljdursi/poapy.get
/*
use std::str;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::collections::HashMap;
use std::collections::HashSet;

use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::banded;
use bio::utils::TextSlice;
use bio::io::fasta;

use petgraph;
use petgraph::Graph;
use petgraph::graph::{NodeIndex, EdgeIndex};
use petgraph::visit::EdgeRef;
use petgraph::dot::Dot;

use util::u8_to_string;
*/
/*
/// A Partial-Order Alignment Graph
#[derive(Default)]
pub struct POAGraph {
    graph: Graph<u8, u16, petgraph::Directed>,
    cns: Vec<u8>,
    cns_path: Vec<NodeIndex>,
    node_idx: Vec<NodeIndex>,
    needs_sort: bool,
}

impl POAGraph {
    /// Create new POAGraph instance, optionally initialized with a sequence
    ///
    /// # Arguments
    ///
    /// * `label` - optional sequence label for an initial sequence
    /// * `sequence` - opetional TextSlice from which to initialize the POA
    ///
    pub fn new_from_sequence(label: Option<&str>, sequence: Option<TextSlice>) -> POAGraph {
        let mut poa = POAGraph {
            graph: Graph::new(),
            cns: Vec::new(),
            cns_path: Vec::new(),
            node_idx: Vec::new(),
            needs_sort: false,
        };

        // Only initialize with a sequence if both a label and sequence were provided
        match (label, sequence) {
            (Some(lab), Some(seq)) => poa.add_unmatched_sequence(lab, seq),
            (_, _) => (None, None),
        };

        return poa;
    }

    /// Create new POAGraph instance, optionally initialized with a sequence
    pub fn new() -> POAGraph {
        POAGraph::new_from_sequence(None::<&str>, None::<TextSlice>)
    }

    /// Return the number of nodes in the underlying graph
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Return the number of edges (note: not edge-weights) in the underlying graph
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Add a new sequence node to the underlying graph
    ///
    /// # Arguments
    ///
    /// * `weight` - The character or nucleotide for this position
    ///
    pub fn add_node(&mut self, weight: u8) -> NodeIndex {
        let new_node = self.graph.add_node(weight);
        self.needs_sort = true;

        new_node
    }

    /// Add a new edge connecting two sequence nodes in the graph
    ///
    /// # Arguments
    ///
    /// * `in_node` - The 'source' of the new edge
    /// * `out_node` - The 'sink' or 'target' of the new edge
    /// * `weight` - The weight associated with the new edge
    ///
    pub fn add_edge(&mut self, in_node: NodeIndex, out_node: NodeIndex, weight: u16) -> EdgeIndex {
        let new_edge = self.graph.add_edge(in_node, out_node, weight);
        self.needs_sort = true;

        new_edge
    }

    /// Add a new sequence new unaligned sequence to the underlying graph.
    /// Useful for both initializing the graph with it's first sequence, as
    /// well as adding unaligned prefix or suffix sequence from partially
    /// aligned reads.
    ///
    /// # Arguments
    ///
    /// * `label` - the id of the sequence to be added
    /// * `sequence` - The sequence to be added
    ///
    pub fn add_unmatched_sequence(
        &mut self,
        _label: &str,
        sequence: TextSlice,
    ) -> (Option<NodeIndex>, Option<NodeIndex>) {
        // Add each character to the graph, and it's associated edges
        let mut first = None::<NodeIndex>;
        let mut last = None::<NodeIndex>;
        for b in sequence {
            let node = self.add_node(*b);
            if first.is_none() {
                first = Some(node);
            }
            if last.is_some() {
                self.add_edge(last.unwrap(), node, 1);
            }
            last = Some(node);
        }

        // Return the book-end positions of the new sequence
        return (first, last);
    }

    /// Calculate and return the current consensus sequence as an array reference
    pub fn consensus(&mut self) -> &[u8] {
        // If we've added no new nodes or edges since the last call, sort first
        if self.needs_sort {

            self.topological_sort();
        }

        self.cns = Vec::new();
        for node in self.consensus_path().to_vec() {
            self.cns.push(*self.graph.node_weight(node).unwrap() as u8);
        }

        &self.cns
    }

    /// Calculate and return the current consensus-path through the underlying
    /// graph structure and return it as an array reference of node indices
    pub fn consensus_path(&mut self) -> &[NodeIndex] {
        // If we've added no new nodes or edges since the last call, sort first
        if self.needs_sort {

            self.topological_sort();
        }

        // For each node find the best predecessor by edge-weight, breaking ties with path-weight
        let mut scores = HashMap::new();
        let mut next_in_path: HashMap<NodeIndex, Option<NodeIndex>> = HashMap::new();
        for node in self.node_idx.iter().rev() {
            let mut best_neighbor = None::<NodeIndex>;
            let mut best_weights = (0, 0); // (Edge-weight, Path-weight)
            for e_ref in self.graph.edges(*node) {
                let weight = *e_ref.weight();
                let target = e_ref.target();

                if (weight, *scores.entry(target).or_insert(0)) > best_weights {
                    best_neighbor = Some(target);
                    best_weights = (weight, *scores.entry(target).or_insert(0));
                }
            }

            scores.insert(*node, best_weights.0 + best_weights.1);
            next_in_path.insert(*node, best_neighbor);
        }

        // Find the start position by locating the highest scoring node
        let mut start_node = None::<NodeIndex>;
        let mut best_score = 0;
        for (node, score) in &scores {
            if score > &best_score {
                start_node = Some(*node);
                best_score = *score;
            }
        }

        // Traverse the graph from the start node, recording the path of nodes taken
        self.cns_path = Vec::new();
        let mut curr_node = start_node;
        //self.write_dot("dot.txt".to_string());
        loop {

            // BUG: Encountering an infinite loop here
            if curr_node.is_some() {
                let node = curr_node.unwrap();
                self.cns_path.push(node);
                curr_node = *next_in_path.get(&node).unwrap();
            } else {
                break;
            }
        }

        &self.cns_path
    }

    /// Align a sequence to the current consensus and return the
    /// resulting alignment object
    ///
    /// # Arguments
    ///
    /// * `sequence` - The new sequence to align as a text slice
    ///
    pub fn align_sequence(&self, seq: &[u8]) -> (Alignment, Vec<NodeIndex>) {

        let match_penalty: i32 = 1;
        let mismatch_penalty: i32 = -1;
        let ins_penalty: i32 = -1;
        let del_penalty: i32 = -1;

        let mut matrix: Vec<Vec<i32>> = vec![vec![0i32; self.node_count()]; seq.len() + 1];
        let mut backtrack: Vec<Vec<Option<(usize, usize)>>> = vec![vec![None; self.node_count()]; seq.len() + 1];

        let mut max_matrix_score = 0;
        let mut max_matrix_loc = (0,0);

        for i in 0..seq.len() {
            for j in 0..self.node_count() {

                let node = self.node_idx[j];
                let node_ix = node.index();
                // take a node. we know it's in topological order
                // so we have definitely already considered any predecessor

                let seq_char: u8 = *seq.get(i).expect("seq[i] shouldn't be empty");

                let node_char: u8 = *self.graph.node_weight(node).unwrap() as u8;

                let match_score = if seq_char == node_char {
                    match_penalty
                } else {
                    mismatch_penalty
                };

                // suppose that seq_char is aligned to node.
                // for each predecessor of node, we need to consider what if we came from that node
                // if we came from that node, then either we're in a match, or an insertion, or a deletion

                let mut max_option= 0;
                let mut max_backtrack: Option<(usize, usize)> = None;
                let mut max_predecessor = 0;

                for incoming in self.graph.neighbors_directed(node, petgraph::Direction::Incoming) {

                    let incoming_ix = incoming.index();
                    // insertion is up
                    // deletion is left
                    // match is diag
                    let ins = matrix[i - 1][node_ix] + ins_penalty;
                    let del = matrix[i][incoming_ix] + del_penalty;
                    let mat = matrix[i - 1][incoming_ix] + match_score;

                    let options: Vec<(i32,(usize,usize))> = vec![(ins,(i-1,node_ix)), (del,(i,incoming_ix)), (mat,(i-1,incoming_ix))];
                    for (option, bt) in options {
                        if option > max_option {
                            max_option = option;
                            max_backtrack = Some(bt);
                        }
                    }
                }

                // local alignment -- allow to start alignment anywhere in the matrix
                if max_option < 0 {
                    max_option = 0;
                    max_backtrack = None;
                }

                // and keep the best scoring spot in the matrix to know where alignment ends
                if max_option > max_matrix_score {
                    max_matrix_score = max_option;
                    max_matrix_loc = (i, node_ix);
                }

                // update the score matrix and backtrack matrix
                matrix[i][node_ix] = max_option;
                backtrack[i][node_ix] = max_backtrack;

            }
        }

        let mut ops: Vec<AlignmentOperation> = vec![];
        // backtrack and store the specific "path" (sequence of DAG nodes) taken
        let mut aln_path: Vec<NodeIndex> = vec![];
        let mut i = max_matrix_loc.0;
        let mut j = max_matrix_loc.1;

        while backtrack[i][j] != None {
            if i - backtrack[i][j].unwrap().0 == 1 && j != backtrack[i][j].unwrap().1 {

                let seq_char: u8 = *seq.get(i).expect("seq[i] shouldn't be empty");
                let node_char: u8 = *self.graph.node_weight(petgraph::prelude::NodeIndex.new(j)).unwrap() as u8;

                if seq_char == node_char {
                    ops.push(AlignmentOperation::Match)
                } else {
                    ops.push(AlignmentOperation::Subst)
                }
            } else if i - backtrack[i][j].unwrap().0 == 1 && j == backtrack[i][j].unwrap().1 {
                ops.push(AlignmentOperation::Ins)
            } else if i == backtrack[i][j].unwrap().0 && j != backtrack[i][j].unwrap().1 {
                ops.push(AlignmentOperation::Del)
            }
            aln_path.push(j);
            (i,j) = backtrack[i][j].unwrap();
        }

        let score = max_matrix_score;
        let ystart = 0;
        let yend = aln_path.len() - 1;
        let xstart = i;
        let xend = max_matrix_loc.0;
        let xlen = seq.len();
        ops.reverse();

        let aln = Alignment {
            score: score,
            ystart: ystart,
            xstart: xstart,
            yend: yend,
            xend: xend,
            xlen: xlen,
            operations: ops,
        };

        // we added the nodes in reverse order so we reverse the path
        aln_path.reverse();

        (aln, aln_path)
    }

    /// Incorporate a new sequence into the graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment object of the new sequence to the current consensus
    /// * `label` - The name of the new sequence being added to the graph
    /// * `seq` - The complete sequence of the read being incorporated
    ///
    pub fn incorporate_alignment(&mut self, aln: Alignment, label: &str, seq: TextSlice) {
        // Ends of sequence may be unaligned, and incorporated separately
        let mut head_node = None::<NodeIndex>;
        let mut tail_node = None::<NodeIndex>;
        if aln.xstart > 0 {
            let (_start_node, end_node) = self.add_unmatched_sequence(&label, &seq[..aln.xstart]);
            head_node = end_node;
        }
        if aln.xend < seq.len() {
            let (start_node, _end_node) = self.add_unmatched_sequence(&label, &seq[aln.xend + 1..]);
            tail_node = start_node;
        }

        // Since each transition is an graph edge, we need to choose an initial source 'prev_node':
        //   1) Use the last node of the unmatched prefix as the initial "prev_node", start at 0
        //   2) Use the node at 0 as the initial "prev_node", start at 1
        let (astart, mut xpos, mut ypos, mut prev_node) = match head_node {
            Some(_n) => (0, aln.xstart, aln.ystart, head_node.unwrap()),
            None => (1, aln.xstart + 1, aln.ystart + 1, self.cns_path[aln.ystart]),
        };

        for apos in astart..aln.operations.len() {
            let op = aln.operations[apos];
            match op {
                AlignmentOperation::Match => {
                    let curr_node = self.cns_path[ypos];

                    match self.graph.find_edge(prev_node, curr_node) {
                        Some(e) => {
                            *self.graph.edge_weight_mut(e).unwrap() += 1;
                        }
                        None => {
                            if self.graph.find_edge(curr_node, prev_node) != None {
                                panic!("Tried to add an edge a->b but edge b->a exists");
                            }
                            self.add_edge(prev_node, curr_node, 1);
                        }
                    }
                    xpos += 1;
                    ypos += 1;
                    prev_node = self.cns_path[ypos - 1];
                }
                AlignmentOperation::Subst => {
                    let xchar = seq[xpos];

                    let out_edge = self.graph
                        .edges(prev_node)
                        .find(|e| *self.graph.node_weight(e.target()).unwrap() as u8 == xchar)
                        .map(|e| e.id());
                    match out_edge {
                        Some(e) => {
                            *self.graph.edge_weight_mut(e).unwrap() += 1;
                            prev_node = self.graph.edge_endpoints(e).unwrap().1;
                            println!("  increasing edge weight");
                        }
                        None => {
                            let new_node = self.add_node(xchar);
                            self.add_edge(prev_node, new_node, 1);
                            prev_node = new_node;
                            println!("  adding new edge");

                        }
                    }
                    xpos += 1;
                    ypos += 1;
                }
                AlignmentOperation::Del => {

                    ypos += 1;
                }
                AlignmentOperation::Ins => {

                    let xchar = seq[xpos];
                    let out_edge = self.graph
                        .edges(prev_node)
                        .filter(|e| e.target() != self.cns_path[ypos])
                        .find(|e| *self.graph.node_weight(e.target()).unwrap() as u8 == xchar)
                        .map(|e| e.id());
                    match out_edge {
                        Some(e) => {
                            *self.graph.edge_weight_mut(e).unwrap() += 1;
                            prev_node = self.graph.edge_endpoints(e).unwrap().1;
                        }
                        None => {
                            let new_node = self.add_node(xchar);
                            self.add_edge(prev_node, new_node, 1);
                            prev_node = new_node;
                        }
                    }
                    xpos += 1;
                }
                _ => ()  // Skip any clipped bases
            }

            // Finally, if the sequence had an unmatched suffix, connect it to the alignment's end
            if tail_node.is_some() {
                self.add_edge(prev_node, tail_node.unwrap(), 1);
            }
        }
    }

    /// Sort the nodes in the underlying graph topologically,
    /// such that every node index preceeds all nodes it has connecting
    /// edges to, and succeeds all nodes that have edges connecting to
    /// it.  This guarantees stable, fast results from traversals
    /// and consensus operations
    pub fn topological_sort(&mut self) {
        let mut sorted_indices: Vec<NodeIndex> = Vec::new();
        let mut completed: HashSet<NodeIndex> = HashSet::new();

        // Depth-first search for unsorted nodes from some defined start-point
        let dfs = |graph: &Graph<u8, u16, petgraph::Directed>,
                   start: NodeIndex,
                   completed: &mut HashSet<NodeIndex>,
                   sorted_indices: &mut Vec<NodeIndex>|
                   -> () {
            let mut stack = vec![start];
            let mut started = HashSet::new();
            loop {
                let node = match stack.pop() {
                    Some(n) => n,
                    None => {
                        break;
                    }
                };

                if completed.contains(&node) {
                    continue;
                }

                if started.contains(&node) {
                    completed.insert(node);
                    sorted_indices.insert(0, node);
                    started.remove(&node);
                    continue;
                }

                let successors = graph
                    .edges(node)
                    .map(|e| e.target())
                    .filter(|n| !completed.contains(n));
                started.insert(node);
                stack.push(node);
                stack.extend(successors);
            }
        };

        // Loop over all uncompleted nodes, performing a depth-first search on each
        loop {

            if sorted_indices.len() == self.node_count() {
                break;
            }

            let mut found = None::<NodeIndex>;
            for node in self.graph.node_indices() {

                if !completed.contains(&node) {
                    found = Some(node);
                    break;
                }
            }

            assert!(found.is_some());

            dfs(
                &self.graph,
                found.unwrap(),
                &mut completed,
                &mut sorted_indices,
            );
        }

        // Finally, store the results of the sort
        self.node_idx = sorted_indices;
        self.needs_sort = false;
    }

    /// Write the current graph to a specified filepath in dot format for
    /// visualization, primarily for debugging / diagnostics
    ///
    /// # Arguments
    ///
    /// * `filename` - The filepath to write the dot file to, as a String
    ///
    pub fn write_dot(&self, filename: String) {
        let mut file = match File::create(&filename) {
            Err(why) => panic!("couldn't open file {}: {}", filename, why.description()),
            Ok(file) => file,
        };
        let g = self.graph.map(|ni, nw| *nw as char,
                               |ni, ew| ew);
        match file.write_all(Dot::new(&g).to_string().as_bytes()) {
            Err(why) => panic!("couldn't write to file {}: {}", filename, why.description()),
            _ => (),
        }
    }
}
*/
pub fn poa_multiple_sequence_alignment(seqlist: &Vec<Vec<u8>>, ref_seq: &[u8]) {

    // obtain reader or fail with error (via ahe unwrap method)

    // Initialize the POA Graph
    let first = ref_seq; //ref_seq.iter().map(|&x| x as u8).collect();
    //let mut poa = POAGraph::new_from_sequence(None, Some(first));
    let mut i = 0;

    for seq in seqlist {
        //let cns = poa.consensus().to_vec().clone();
        //println!("As vec {:?}", cns);
        //println!("As str {}", String::from_utf8(cns.to_vec()).unwrap());
        /*println!(
            "Graph has {} nodes and {} edges",
            poa.node_count(),
            poa.edge_count()
        );*/

        // obtain record or fail with error
        //println!("\n{}: {}", i, record.id());

        //let alignment = poa.align_sequence(&seq);
        //println!("{}", alignment.pretty(record.seq(), &cns));
        println!(">seq{}",i);
        println!("{}",u8_to_string(&seq));
        //poa.incorporate_alignment(alignment, &i.to_string(), &seq);
        i += 1
    }

    //u8_to_string(poa.consensus().clone())

    //println!("As vec {:?}", cns);
    //println!("As str {}", String::from_utf8(cns.to_vec()).unwrap());
    //println!(
    //    "Graph has {} nodes and {} edges",
    //    poa.node_count(),
    //    poa.edge_count()
    //);

    //poa.write_dot("test.dot".to_string());

    //println!("{}",));
    //cns
}
/*
pub fn test() {

    let filename = "poa_failure_case_cycle.fa";

    // obtain reader or fail with error (via ahe unwrap method)

    let file = File::open(filename).unwrap();
    let reader = fasta::Reader::new(file);
    let mut records = reader.records();

    // Initialize the POA Graph
    let first = records.next().unwrap().unwrap();
    let mut poa = POAGraph::new_from_sequence(Some(first.id()), Some(first.seq()));

    for (i, result) in records.enumerate() {
        let cns = poa.consensus().to_vec().clone();
        println!("As vec {:?}", cns);
        println!("As str {}", String::from_utf8(cns.to_vec()).unwrap());
        println!(
            "Graph has {} nodes and {} edges",
            poa.node_count(),
            poa.edge_count()
        );

        // obtain record or fail with error
        let record = result.unwrap();
        println!("\n{}: {}", i, record.id());

        let alignment = poa.align_sequence(record.seq());
        println!("{}", alignment.pretty(record.seq(), &cns));
        poa.incorporate_alignment(alignment, record.id(), record.seq());
    }

    let cns = poa.consensus().to_vec().clone();
    println!("As vec {:?}", cns);
    println!("As str {}", String::from_utf8(cns.to_vec()).unwrap());
    println!(
        "Graph has {} nodes and {} edges",
        poa.node_count(),
        poa.edge_count()
    );

    poa.write_dot("test.dot".to_string());
}
*/
