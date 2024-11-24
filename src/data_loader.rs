// Modules for data loading
use std::collections::{HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

// Function to load the cell barcodes
pub fn load_cell_barcodes(file_path: Option<&String>) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    let mut barcodes = HashSet::new();
    if let Some(path) = file_path {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let barcode = line?.trim().to_string();
            barcodes.insert(barcode);
        }
    }
    Ok(barcodes)
}
