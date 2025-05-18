use std::{cmp::Ordering, slice};

fn main() {
    // Our reference genome, stored in a Vec<u8>
    let reference = b"ATGCGT".to_vec();

    // Build a naive suffix array by comparing substrings via raw pointers
    let suffix_array = build_suffix_array(&reference);
    println!("Suffix Array for {:?}: {:?}", String::from_utf8_lossy(&reference), suffix_array);
}

// Given a vector of bytes, we generate a suffix array by sorting indices [0..len]
// according to the lexicographical order of the suffixes they represent.
fn build_suffix_array(seq: &Vec<u8>) -> Vec<usize> {
    let len = seq.len();
    let ptr = seq.as_ptr(); // raw pointer to the first byte of seq

    // We collect all starting positions [0..len].
    let mut suffixes = (0..len).collect::<Vec<usize>>();

    // Sort the indices by comparing suffixes using a custom comparator.
    // We must resort to 'unsafe' to manipulate raw pointers,
    // but Rust ensures that these pointers are valid within seq's lifetime.
    suffixes.sort_by(|&i, &j| unsafe {
        compare_suffixes(ptr, len, i, j)
    });

    suffixes
}

// Compare two suffixes using raw pointers and lexicographical ordering.
unsafe fn compare_suffixes(ptr: *const u8, len: usize, i: usize, j: usize) -> Ordering {
    // Construct slices from raw pointers, ensuring we only compare valid memory.
    let suffix_i = slice::from_raw_parts(ptr.add(i), len - i);
    let suffix_j = slice::from_raw_parts(ptr.add(j), len - j);

    suffix_i.cmp(suffix_j)
}