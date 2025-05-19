use rust_htslib::bcf;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

pub fn phase_block(record: &mut bcf::Record, _window_size: i32) -> Result<String, String> {
    // This is a mock implementation that just assigns a phase block based on chromosome and position
    // Get the chromosome name (using rid and header instead of chrom())
    let rid = record.rid().ok_or_else(|| "Record has no RID".to_string())?;
    let header = record.header();
    let chrom = std::str::from_utf8(header.rid2name(rid).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    
    // Get position
    let pos = record.pos();
    
    // Create a deterministic phase block ID based on chromosome and position range
    let phase_block_id = format!("{}:{}", chrom, pos / 1000 * 1000);
    
    // Hash the phase block ID to create a shorter identifier
    let mut hasher = DefaultHasher::new();
    phase_block_id.hash(&mut hasher);
    let hash_value = hasher.finish();
    
    Ok(format!("PB{:X}", hash_value % 1000))
}
