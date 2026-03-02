#!/usr/bin/env Rscript

library(alakazam)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: Rscript hill_diversity_single.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

# Ensure output directory exists
out_dir <- dirname(output_file)
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

message("Processing ", input_file)

# Load CSV
data <- read.csv(input_file)

# Run Hill diversity profile
hill <- alphaDiversity(
    data = data,
    group = "patient_id",
    clone = "counts",
    min_q = 0,
    max_q = 10,
    step_q = 0.1,
    ci = 0.95,
    nboot = 2000
)

# Save output
write_tsv(hill@diversity, output_file)

message("Saved: ", output_file)
