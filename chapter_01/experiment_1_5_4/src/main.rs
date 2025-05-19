use polars::prelude::*;
use std::error::Error;

fn read_csv_basic() -> Result<DataFrame, Box<dyn Error>> {
    let df_csv = LazyCsvReader::new("output.csv")
        .with_has_header(true) // ✅ FIXED
        .finish()?
        .collect()?; // Convert LazyFrame to DataFrame
    
    Ok(df_csv)
}

fn read_csv_advanced() -> Result<DataFrame, Box<dyn Error>> {
    let df_csv = LazyCsvReader::new("output.csv")
        .with_has_header(true) // ✅ FIXED
        .with_separator(b',') // ✅ FIXED: Replaced `with_delimiter()` with `with_separator()`
        .with_null_values(Some(NullValues::AllColumns(vec!["NA".into(), "null".into()])))
        .finish()?
        .collect()?; // Convert LazyFrame to DataFrame
    
    Ok(df_csv)
}

fn main() -> Result<(), Box<dyn Error>> {
    // Try reading with basic options
    let df = read_csv_basic()?;
    println!("Basic CSV read result:");
    println!("{}", df.head(Some(5))); // Print first 5 rows

    // Try reading with advanced options
    let df_advanced = read_csv_advanced()?;
    println!("\nAdvanced CSV read result:");
    println!("{}", df_advanced.head(Some(5))); // Print first 5 rows

    // Ensure DataFrame has data before grouping
    if df.height() > 0 {
        let lazy_df = df.lazy()
            .group_by([col("disease_state")])
            .agg([
                col("quality_score").mean().alias("avg_quality"),
                col("patient_age").mean().alias("avg_age"),
                col("sample_id").count().alias("sample_count"),
            ])
            .collect()?; // Convert LazyFrame back to DataFrame

        println!("\nSummary by disease state:");
        println!("{}", lazy_df);
    }

    Ok(())
}
