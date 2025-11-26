#!/usr/bin/env Rscript

## LAVA local rg analysis (univariate + bivariate) for a block of loci
## Usage:
## Rscript run_lava_local_rg.R \
##   <ref_prefix> <locfile> <input.info> <sample.overlap> <phenos> <out_prefix> <start> <stop>
##
## Example:
## Rscript run_lava_local_rg.R \
##   data/refdat/g1000_eur.maf01.ids \
##   data/loci/blocks_s2500_m25_f1_w200.locfile \
##   data/input.info.txt \
##   data/sample.overlap.txt \
##   DCM;HCM \
##   results/lava_dcm_hcm_out_ \
##   1 5

suppressPackageStartupMessages({
  library(parallel)
  library(LAVA)
  library(snpStats)
  library(data.table)
  library(matrixsampling)
  library(keep)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  stop(
    "Expected 8 arguments:\n",
    " 1) ref_prefix (e.g. data/refdat/g1000_eur.maf01.ids)\n",
    " 2) locfile (e.g. data/loci/blocks_s2500_m25_f1_w200.locfile)\n",
    " 3) input.info (e.g. data/input.info.txt)\n",
    " 4) sample.overlap (e.g. data/sample.overlap.txt)\n",
    " 5) phenos (e.g. DCM;HCM)\n",
    " 6) out_prefix (e.g. results/lava_dcm_hcm_out_)\n",
    " 7) start_locus (integer)\n",
    " 8) stop_locus (integer)\n",
    call. = FALSE
  )
}

ref_prefix   <- args[1]
locfile      <- args[2]
input_info   <- args[3]
overlap_file <- args[4]
phenos       <- unlist(strsplit(args[5], ";"))
out_prefix   <- args[6]
start_locus  <- as.integer(args[7])
stop_locus   <- as.integer(args[8])

if (any(is.na(c(start_locus, stop_locus)))) {
  stop("start_locus and stop_locus must be integers.", call. = FALSE)
}

univ_thresh <- 0.05 / 2495  # example: 0.05 / N_loci (adjust as needed)

## ---- Read loci definition ----
loci <- LAVA::read.loci(locfile)
n_loci_total <- nrow(loci)

if (stop_locus > n_loci_total) {
  stop("stop_locus exceeds number of loci in locfile: ", n_loci_total, call. = FALSE)
}

message("Running LAVA for loci ", start_locus, " to ", stop_locus, " out of ", n_loci_total)

## ---- Process input (GWAS + LD + sample overlap) ----
input <- LAVA::process.input(
  input.info.file     = input_info,
  sample.overlap.file = overlap_file,
  ref.prefix          = ref_prefix,
  phenos              = phenos
)

## ---- Analysis loop ----
loc_indices <- start_locus:stop_locus
progress    <- ceiling(stats::quantile(loc_indices, seq(0.05, 1, 0.05)))

univ_results  <- NULL
bivar_results <- NULL

message(Sys.time(), "  Starting LAVA analysis for ", length(loc_indices), " loci")

for (i in loc_indices) {

  if (i %in% progress) {
    pct <- names(progress[which(progress == i)])
    message(Sys.time(), "  .. ", pct, " of loci completed")
  }

  locus_obj <- LAVA::process.locus(loci[i, ], input)

  # locus may be NULL (too few SNPs, negative variances, etc.)
  if (!is.null(locus_obj)) {

    loc_info <- data.frame(
      locus  = locus_obj$id,
      chr    = locus_obj$chr,
      start  = locus_obj$start,
      stop   = locus_obj$stop,
      n.snps = locus_obj$n.snps,
      n.pcs  = locus_obj$K
    )

    ub <- LAVA::run.univ.bivar(locus_obj, univ.thresh = univ_thresh)

    if (!is.null(ub$univ)) {
      univ_results <- rbind(univ_results, cbind(loc_info, ub$univ))
    }

    if (!is.null(ub$bivar)) {
      bivar_results <- rbind(bivar_results, cbind(loc_info, ub$bivar))
    }
  }
}

## ---- Write output ----
if (!dir.exists(dirname(out_prefix))) {
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
}

univ_file  <- paste0(out_prefix, start_locus, "_", stop_locus, ".univ")
bivar_file <- paste0(out_prefix, start_locus, "_", stop_locus, ".univFDR5.bivar")

if (!is.null(univ_results)) {
  write.table(univ_results, univ_file,
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  message("Univariate results written to: ", univ_file)
} else {
  message("No univariate results produced for this locus range.")
}

if (!is.null(bivar_results)) {
  write.table(bivar_results, bivar_file,
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  message("Bivariate results written to: ", bivar_file)
} else {
  message("No bivariate results produced for this locus range.")
}

message(Sys.time(), "  LAVA analysis finished for loci ", start_locus, "-", stop_locus)
