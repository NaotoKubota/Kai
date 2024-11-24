[![GitHub License](https://img.shields.io/github/license/NaotoKubota/Kai)](https://github.com/NaotoKubota/Kai/blob/main/LICENSE)
[![GitHub Release](https://img.shields.io/github/v/release/NaotoKubota/Kai?style=flat)](https://github.com/NaotoKubota/Kai/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/NaotoKubota/Kai)](https://github.com/NaotoKubota/Kai/releases)
[![Rust](https://github.com/NaotoKubota/Kai/actions/workflows/rust.yaml/badge.svg)](https://github.com/NaotoKubota/Kai/actions/workflows/rust.yaml)
[![Create Release and Build Docker Image](https://github.com/NaotoKubota/Kai/actions/workflows/release-docker-build-push.yaml/badge.svg)](https://github.com/NaotoKubota/Kai/actions/workflows/release-docker-build-push.yaml)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/kai)](https://hub.docker.com/r/naotokubota/kai)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/kai)](https://hub.docker.com/r/naotokubota/kai)

# Kai (v0.1.0)

Rust implementation of read counting for regions of interest from RNA-seq BAM files.

## Usage

```bash
Usage: kai [OPTIONS] <mode> <bam_file> <regions_file> <output_prefix>

Arguments:
  <mode>           Mode of operation: 'bulk' or 'single' [possible values: bulk, single]
  <bam_file>       Path to the BAM file
  <regions_file>   Path to the BED file containing regions of interest
  <output_prefix>  Output prefix for the output files

Options:
  -l, --max-loci <max_loci>                Maximum number of loci the read maps to [default: 1]
  -c, --cell-barcodes <cell_barcode_file>  Optional file specifying cell barcodes of interest
  -v, --verbose                            Enable verbose output to print all arguments
  -h, --help                               Print help
  -V, --version                            Print version
```

## Build

```bash
cargo build --release
```

## Example

```bash
# Count junction reads from bulk RNA-seq BAM file
./target/release/kai bulk example.bam regions.bed output_example
# Count junction reads from single-cell RNA-seq BAM file
./target/release/kai single example.bam regions.bed output_example
```
