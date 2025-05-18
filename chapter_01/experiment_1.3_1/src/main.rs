fn build_prefix_table(pattern: &str) -> Vec<usize> {
    let m = pattern.len();
    let mut pi = vec![0; m];
    let mut j = 0;
    let bytes = pattern.as_bytes();

    for i in 1..m {
        while j > 0 && bytes[i] != bytes[j] {
            j = pi[j - 1];
        }
        if bytes[i] == bytes[j] {
            j += 1;
            pi[i] = j;
        }
    }
    pi
}

fn kmp_search(text: &str, pattern: &str) -> Vec<usize> {
    let n = text.len();
    let m = pattern.len();
    if m == 0 || n < m {  // Cek jika pola kosong atau lebih panjang dari teks
        return vec![];
    }
    let pi = build_prefix_table(pattern);
    let mut matches = Vec::new();

    let t_bytes = text.as_bytes();
    let p_bytes = pattern.as_bytes();
    let mut j = 0;

    for i in 0..n {
        while j > 0 && t_bytes[i] != p_bytes[j] {
            j = pi[j - 1];
        }
        if t_bytes[i] == p_bytes[j] {
            j += 1;
        }
        if j == m {
            if i >= m - 1 { // Cek untuk menghindari underflow
                matches.push(i + 1 - m);
            }
            j = pi[j - 1];
        }
    }
    matches
}

fn main() {
    let dna_sequence = "ATGCGATATCGATGCGATGCGATGC";
    let pattern = "ATGC";

    let matches = kmp_search(dna_sequence, pattern);

    if matches.is_empty() {
        println!("Pola '{}' tidak ditemukan dalam DNA sequence.", pattern);
    } else {
        println!("Pola '{}' ditemukan pada indeks: {:?}", pattern, matches);
    }
}
