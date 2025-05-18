pub mod graph {
    use serde::{Deserialize, Serialize};
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;

    #[derive(Debug, Serialize, Deserialize)]
    struct NodeData {
        id: u64,
        sequence: String,
        chrom: String,
        pos: u64,
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct EdgeData {
        from: u64,
        to: u64,
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct GraphData {
        nodes: Vec<NodeData>,
        edges: Vec<EdgeData>,
        metadata: HashMap<String, String>,
    }

    #[derive(Debug)]
    pub struct Graph {
        nodes: HashMap<u64, NodeData>,
        edges: Vec<EdgeData>,
        chrom_pos_map: HashMap<(String, u64), u64>, // (chrom, pos) -> node_id
    }

    impl Graph {
        pub fn from_json_path<P: AsRef<Path>>(path: P) -> Result<Self, String> {
            let file = File::open(path).map_err(|e| e.to_string())?;
            let reader = BufReader::new(file);
            
            let data: GraphData = serde_json::from_reader(reader)
                .map_err(|e| format!("Failed to parse JSON: {}", e))?;
            
            let mut nodes = HashMap::new();
            let mut chrom_pos_map = HashMap::new();
            
            for node in data.nodes {
                chrom_pos_map.insert((node.chrom.clone(), node.pos), node.id);
                nodes.insert(node.id, node);
            }
            
            Ok(Graph {
                nodes,
                edges: data.edges,
                chrom_pos_map,
            })
        }
        
        pub fn node_count(&self) -> usize {
            self.nodes.len()
        }
        
        pub fn edge_count(&self) -> usize {
            self.edges.len()
        }
        
        pub fn node_at(&self, chrom: &str, pos: u64) -> Option<u64> {
            self.chrom_pos_map.get(&(chrom.to_string(), pos)).copied()
        }
        
        pub fn degree_at(&self, chrom: &str, pos: u64) -> Option<u32> {
            self.node_at(chrom, pos).map(|id| self.degree(id))
        }
        
        pub fn centrality_at(&self, chrom: &str, pos: u64) -> Option<f64> {
            self.node_at(chrom, pos).map(|id| self.centrality(id))
        }
        
        pub fn degree(&self, node_id: u64) -> u32 {
            let mut count = 0;
            for edge in &self.edges {
                if edge.from == node_id || edge.to == node_id {
                    count += 1;
                }
            }
            count
        }
        
        pub fn centrality(&self, _node_id: u64) -> f64 {
            // Mock implementation: Just return a random value between 0 and 1
            rand::random::<f64>()
        }
    }
}
