use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use rayon::prelude::*;
use bio::io::fasta;

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct PGNode {
    kmer: String,
}

type PGGraph = HashMap<PGNode, Vec<PGNode>>;

fn build_pangenome_graph(k: usize, fasta_files: &[&str]) -> PGGraph {
    let partial_graphs: Vec<PGGraph> = fasta_files
        .par_iter()
        .map(|path| {
            let reader = fasta::Reader::from_file(path)
                .unwrap_or_else(|_| panic!("Cannot open FASTA file: {}", path));
            
            let mut local_graph = PGGraph::new();
            for record in reader.records() {
                let rec = record.expect("Invalid FASTA record");
                let seq_bytes = rec.seq();
                if seq_bytes.len() < k {
                    continue;
                }
                let seq_str = String::from_utf8(seq_bytes.to_vec())
                    .expect("Non-UTF8 sequence data");
                
                for i in 0..seq_str.len().saturating_sub(k) {
                    let node_kmer = &seq_str[i..i + k];
                    let edge_kmer = &seq_str[i + 1..i + k + 1];

                    let node = PGNode {
                        kmer: node_kmer.to_owned(),
                    };
                    let edge_node = PGNode {
                        kmer: edge_kmer.to_owned(),
                    };

                    local_graph
                        .entry(node)
                        .or_insert_with(Vec::new)
                        .push(edge_node);
                }
            }
            local_graph
        })
        .collect();

    partial_graphs.into_par_iter().reduce(
        || PGGraph::new(),
        |mut acc, local| {
            for (node, successors) in local {
                acc.entry(node)
                    .or_insert_with(Vec::new)
                    .extend(successors);
            }
            acc
        },
    )
}

fn main() {
    let haplotypes = &["src/haplotype1.fasta", "src/haplotype2.fasta", "src/haplotype3.fasta"];
    let k = 21;
    let pangenome_graph = build_pangenome_graph(k, haplotypes);

    let mut file = File::create("output.txt").expect("Cannot create output file");
    writeln!(file, "Constructed a pangenome graph with {} nodes.", pangenome_graph.len())
        .expect("Cannot write to file");
    
    for (node, edges) in pangenome_graph.iter().take(5) {
        writeln!(file, "Node: {} -> {:?}", node.kmer, edges.iter().map(|e| &e.kmer).collect::<Vec<_>>())
            .expect("Cannot write to file");
    }
}
