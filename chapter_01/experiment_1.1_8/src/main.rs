fn main() {
    // A collection of short DNA reads for demonstration purposes
    let reads = vec!["ATCG", "TTGA", "CGTA", "GGTT"];

    // We transform each read (map) by counting GC bases, 
    // then combine (fold) these counts into a single total.
    let total_gc = reads
        .iter()
        .map(|read| {
            read.chars()
                .filter(|&c| c == 'G' || c == 'C')
                .count()
        })
        .fold(0, |acc, gc_count| acc + gc_count);

    println!("Total GC content in all reads: {}", total_gc);
}
