use mpi::traits::*;
use mpi::datatype::PartitionMut;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;

/// Simulates a local index chunk, for example, a partial FM-index or suffix array for a slice of the genome.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct PartialIndex {
    rank_id: i32,
    index_data: HashMap<String, usize>,
}

/// Example function to build an index for a slice of the genome on one HPC node.
fn build_local_index(rank: i32) -> PartialIndex {
    let mut data = HashMap::new();
    data.insert(format!("key_rank_{}", rank), rank as usize);
    PartialIndex {
        rank_id: rank,
        index_data: data,
    }
}

/// Merges all partial indexes into a single global index.
fn merge_indexes(chunks: Vec<PartialIndex>) -> HashMap<String, usize> {
    let mut global_map = HashMap::new();
    for chunk in chunks {
        for (k, v) in chunk.index_data {
            global_map.insert(k, v);
        }
    }
    global_map
}

fn main() {
    // Initialize the MPI environment
    let universe = mpi::initialize().expect("Failed to initialize MPI environment");
    let world = universe.world();
    let rank = world.rank();
    let size = world.size();

    // Each rank builds a local partial index
    let local_chunk = build_local_index(rank);

    // Serialize the data for transmission
    let serialized_local_chunk = bincode::serialize(&local_chunk).expect("Serialization failed");

    // Gather the lengths of the serialized chunks from all ranks
    let local_len = serialized_local_chunk.len();
    let mut recv_lengths = vec![0; size as usize];
    world.all_gather_into(&local_len, &mut recv_lengths);

    // Convert recv_lengths to i32 slices (fixing the type issue)
    let recv_lengths_i32: Vec<i32> = recv_lengths.iter().map(|&x| x as i32).collect();
    let recv_lengths_slice: &[i32] = &recv_lengths_i32; // Convert Vec<i32> to slice

    // Compute total received buffer size
    let total_recv_size: usize = recv_lengths.iter().sum();
    let mut recv_data = vec![0u8; total_recv_size];

    // Compute displacements (where each rank's data starts)
    let mut displacements = vec![0; size as usize];
    for i in 1..size as usize {
        displacements[i] = displacements[i - 1] + recv_lengths[i - 1];
    }

    // Convert displacements to i32 slices (fixing the type issue)
    let displacements_i32: Vec<i32> = displacements.iter().map(|&x| x as i32).collect();
    let displacements_slice: &[i32] = &displacements_i32; // Convert Vec<i32> to slice

    // Use `PartitionMut` for correct MPI communication
    let mut recv_partition = PartitionMut::new(&mut recv_data, recv_lengths_slice, displacements_slice);
    world.all_gather_varcount_into(&serialized_local_chunk[..], &mut recv_partition);

    // Rank 0 processes the gathered data
    if rank == 0 {
        let mut gathered = Vec::new();
        let mut offset = 0;

        for &len in &recv_lengths {
            let bytes = &recv_data[offset..offset + len];
            let chunk: PartialIndex = bincode::deserialize(bytes).expect("Deserialization failed");
            gathered.push(chunk);
            offset += len;
        }

        // Merge all partial indexes
        let global_index = merge_indexes(gathered);
        println!("Final merged index on rank 0: {:?}", global_index);
    }
}
