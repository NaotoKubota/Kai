use clap::{Arg, Command};
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::bam::record::{Aux, Cigar};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Write, BufRead};
use std::io;
use log::{info, debug, LevelFilter};
use env_logger;
use itertools::Itertools;
use flate2::write::GzEncoder;
use flate2::Compression;

mod data_loader;

struct Region {
    chromosome: String,
    start: usize,
    end: usize,
}

fn parse_bed_file(bed_file: &str) -> io::Result<Vec<Region>> {
    let file = File::open(bed_file)?;
    let reader = io::BufReader::new(file);

    let mut regions = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue; // Skip headers or empty lines
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            debug!("Invalid BED format: {}", line);
            continue;
        }

        let chromosome = fields[0].to_string();
        let start = fields[1].parse::<usize>().expect("Invalid start coordinate");
        let end = fields[2].parse::<usize>().expect("Invalid end coordinate");

        regions.push(Region {
            chromosome,
            start,
            end,
        });
    }

    info!("Parsed {} regions", regions.len());
    Ok(regions)

}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Set up command-line arguments using clap
    let matches = Command::new("kai")
        .version("0.1.1")
        .author("NaotoKubota")
        .about("Count reads mapped to regions of interest from bulk/single-cell RNA-seq data")
        .arg(Arg::new("mode")
            .required(true)
            .value_parser(["bulk", "single"])
            .help("Mode of operation: 'bulk' or 'single'"))
        .arg(Arg::new("bam_file")
            .required(true)
            .help("Path to the BAM file"))
        .arg(Arg::new("regions_file")
            .required(true)
            .help("Path to the BED file containing regions of interest"))
        .arg(Arg::new("output_prefix")
            .required(true)
            .help("Output prefix for the output files"))
        .arg(Arg::new("max_loci")
            .short('l')
            .long("max-loci")
            .default_value("1")
            .value_parser(clap::value_parser!(u32))
            .help("Maximum number of loci the read maps to"))
        .arg(Arg::new("cell_barcode_file")
            .short('c')
            .long("cell-barcodes")
            .value_parser(clap::value_parser!(String))
            .help("Optional file specifying cell barcodes of interest"))
        .arg(Arg::new("verbose")
            .short('v')
            .long("verbose")
            .action(clap::ArgAction::SetTrue)
            .help("Enable verbose output to print all arguments"))
        .get_matches();

    // Parse arguments
    let mode = matches.get_one::<String>("mode").unwrap();
    let bam_file = matches.get_one::<String>("bam_file").unwrap();
    let regions_file = matches.get_one::<String>("regions_file").unwrap();
    let output_prefix = matches.get_one::<String>("output_prefix").unwrap();
    let max_loci = *matches.get_one::<u32>("max_loci").unwrap();
    let cell_barcode_file = matches.get_one::<String>("cell_barcode_file");
    let verbose = matches.get_flag("verbose");

    // Initialize the logger with the appropriate level
    if verbose {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Debug)
            .init();
    } else {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Info)
            .init();
    }

    // Log all arguments if verbose is enabled
    info!("Running kai");
    info!("Mode: {}", mode);
    info!("BAM file: {}", bam_file);
    info!("Regions file: {}", regions_file);
    info!("Output prefix: {}", output_prefix);
    info!("Maximum loci (NH): {}", max_loci);
    // Load cell barcodes of interest
    let cell_barcodes_of_interest = if mode == "single" {
        let barcodes = data_loader::load_cell_barcodes(cell_barcode_file)?;
        info!(
            "Cell barcodes of interest: {}",
            if barcodes.is_empty() {
                "None (processing all reads)".to_string()
            } else {
                format!("{} barcodes", barcodes.len())
            }
        );
        barcodes
    } else {
        HashSet::new()
    };

    // Parse the BED file containing regions of interest
    info!("Parsing regions of interest from BED file");
    let regions = parse_bed_file(&regions_file)?;

    // Prepare a map for counting reads per region and optionally by cell barcode
    let mut region_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut region_totals: HashMap<String, u32> = HashMap::new();
    let mut cell_barcodes: HashSet<String> = HashSet::new();

    // Open the BAM index
    let mut bam = IndexedReader::from_path(bam_file)?;

    // Counter for tracking the number of regions processed
    let mut region_counter = 0;
    let mut last_percentage = 0;

    // Count reads mapped to regions of interest
    info!("Counting reads mapped to regions of interest");
    for region in &regions {
        let region_key = format!("{}:{}-{}", region.chromosome, region.start, region.end);
        region_counter += 1;

        // Calculate and log progress at each 1% increment
        let progress_percentage = (region_counter * 100) / regions.len();
        if progress_percentage > last_percentage {
            info!("Progress: {}% / ({} / {})", progress_percentage, region_counter, regions.len());
            last_percentage = progress_percentage;
        }

        // Fetch reads in the region
        let chrom_bytes = region.chromosome.as_bytes();
        bam.fetch((chrom_bytes, region.start as i64, region.end as i64))?;

        // Iterate over reads in the region
        for result in bam.records() {
            let record = result?;
            // Skip read if NH tag exceeds max_loci
            if let Ok(Aux::U8(nh)) = record.aux(b"NH") {
                if nh > max_loci as u8 {
                    continue; // Skip reads with more than max_loci loci
                }
            }

            // Extract Cell Barcode (CB) from tags if in single mode
            let cell_barcode = if mode == "single" {
                match record.aux(b"CB") {
                    Ok(Aux::String(cb)) => Some(cb.to_string()),
                    _ => None,
                }
            } else {
                None
            };

            // Skip read if its barcode is not in the list of interest
            if let Some(cb) = &cell_barcode {
                if cell_barcode_file.is_some() && !cell_barcodes_of_interest.is_empty() && !cell_barcodes_of_interest.contains(cb) {
                    continue; // Skip reads with cell barcodes not in the list of interest
                }
                cell_barcodes.insert(cb.clone());
            }

            // Get the start position of the read
            let mut current_pos = record.pos();

            // Check cigar string to determine if the read overlaps the region with matching bases, not like RefSkip or SoftClip
            let cigar_vec = record.cigar(); // Create a longer-lived binding for the cigar data
            let cigars: Vec<_> = cigar_vec.iter().collect();
            for i in 0..cigars.len() {
                let cigar = cigars[i];
                if let Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) = cigar {
                    let cigar_len = cigar.len() as i64;
                    let cigar_end = current_pos + cigar_len;
                    // Increment the count if the read overlaps the region at least partially
                    if current_pos < region.end.try_into().unwrap() && cigar_end > region.start.try_into().unwrap() {
                        if mode == "single" {
                            if let Some(cb) = &cell_barcode {
                                let region_entry = region_counts
                                    .entry(region_key.to_string())
                                    .or_insert_with(HashMap::new);
                                *region_entry.entry(cb.clone()).or_insert(0) += 1;
                            }
                        } else if mode == "bulk" {
                            *region_totals
                            .entry(region_key.to_string())
                            .or_insert(0) += 1;
                        }
                        break; // Break the loop to avoid double counting
                    }
                } else if let Cigar::SoftClip(_) = cigar {
                    continue;
                } else {
                    current_pos += match cigar {
                        Cigar::Ins(l) | Cigar::Del(l) | Cigar::RefSkip(l) => *l as i64,
                        _ => 0,
                    };
                }
            }
        }
    }

    // Write results based on mode
    info!("Writing output files");
    if mode == "single" {
        // Prepare output files with compression
        let mut matrix_file = GzEncoder::new(File::create(format!("{}_matrix.mtx.gz", output_prefix))?, Compression::default());
        let mut barcodes_file = GzEncoder::new(File::create(format!("{}_barcodes.tsv.gz", output_prefix))?, Compression::default());
        let mut features_file = GzEncoder::new(File::create(format!("{}_features.tsv.gz", output_prefix))?, Compression::default());
        let mut output_tsv = GzEncoder::new(File::create(format!("{}_count_barcodes.tsv.gz", output_prefix))?, Compression::default());

        // Write barcodes.tsv.gz
        debug!("Writing barcodes.tsv.gz");
        let barcode_list: Vec<_> = cell_barcodes.iter().sorted().collect();
        for barcode in &barcode_list {
            writeln!(barcodes_file, "{}", barcode)?;
        }

        // Write features.tsv.gz
        debug!("Writing features.tsv.gz");
        let feature_list: Vec<_> = region_counts.keys().sorted().collect();
        for feature in &feature_list {
            writeln!(features_file, "{}", feature)?;
        }

        // Buffers to accumulate lines for matrix.mtx.gz and output.tsv.gz
        let mut matrix_buffer: Vec<String> = Vec::new();
        let mut tsv_buffer: Vec<String> = Vec::new();

        // Add the header lines to the matrix buffer
        matrix_buffer.push("%%MatrixMarket matrix coordinate integer general".to_string());
        matrix_buffer.push("%".to_string());
        matrix_buffer.push(format!(
            "{} {} {}",
            feature_list.len(),
            barcode_list.len(),
            region_counts.values().map(|c| c.len()).sum::<usize>()
        ));

        // Add sparse matrix data and TSV data to the buffers
        debug!("Writing matrix.mtx.gz and count_barcodes.tsv.gz");
        let barcode_map: HashMap<_, _> = barcode_list.iter().enumerate().map(|(i, b)| (b.as_str(), i + 1)).collect();
        tsv_buffer.push("Feature\tBarcode\tCount".to_string());
        for (i, feature) in feature_list.iter().enumerate() {
            if let Some(cell_counts) = region_counts.get(*feature) {
                for (barcode, count) in cell_counts {
                    if let Some(&j) = barcode_map.get(barcode.as_str()) {
                        matrix_buffer.push(format!("{} {} {}", i + 1, j + 1, count));
                        tsv_buffer.push(format!("{}\t{}\t{}", feature, barcode, count));
                    }
                }
            }
        }

        // Write the accumulated lines to the compressed output files
        for line in matrix_buffer {
            writeln!(matrix_file, "{}", line)?;
        }
        for line in tsv_buffer {
            writeln!(output_tsv, "{}", line)?;
        }

    } else {
        let mut output_file = GzEncoder::new(File::create(format!("{}_count.tsv.gz", output_prefix))?, Compression::default());
        debug!("Writing count.tsv.gz");
        writeln!(output_file, "Chr\tStart\tEnd\tRegion\tCount")?;
        for (region, count) in region_totals.iter().sorted() {
            let fields: Vec<&str> = region.split([':', '-']).collect();
            writeln!(output_file, "{}\t{}\t{}\t{}", fields[0], fields[1], fields[2], count)?;
        }
    }

    info!("Finished processing");
    Ok(())
}
