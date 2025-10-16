options(scipen = 999)
suppressPackageStartupMessages(library(optparse))  # Parsing command-line options
suppressPackageStartupMessages(library(data.table))  # Efficient data reading and manipulation
suppressPackageStartupMessages(library(dplyr))  # Data manipulation and transformation
suppressPackageStartupMessages(library(readr))  # Reading and writing data
suppressPackageStartupMessages(library(tidyr))  # Data tidying

print(paste("Libraries loaded successfully at", Sys.time()))

# Read the TSV file
option_list <- list(
  make_option(c("-f", "--file"), type="character", default="data/raw/Garnier_Meder_Amsterdam_FinnGen_UKB_MGB__DCM__META1_chr1_22_MAF0.005.tsv", 
              help="Path to the input sumstats file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="data/sumst_processed/DCM_GWAS_37_exclMYBPC3reg.txt",
              help="Path to the output file", metavar="file")
)

print("Command-line options set up")

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))
print(paste("Input file:", opt$file))
print(paste("Output file:", opt$output))

# Read the TSV file
print("Reading data from input file")
dat1 <- fread(opt$file)
dat1 <- fread("data/raw/Garnier_Meder_Amsterdam_FinnGen_UKB_MGB__DCM__META1_chr1_22_MAF0.005.tsv")

print(paste("Data read successfully. Rows:", nrow(dat1), ", Columns:", ncol(dat1)))


# Process the data: select and rename columns
print("Reading data from input file")
dat <- dat1

print("Processing data: calculating minor allele frequency (maf)")
dat$maf <- dat$EAFREQ

dat[dat$EAFREQ > 0.5, 'maf'] <- 1 - dat[dat$EAFREQ > 0.5, 'maf']

print(paste("EAFREQ > 0.5. Rows:", nrow(dat), ", Columns:", ncol(dat)))

print("Calculating effective sample size (Neff)")
dat$Neff <- 1 / (2 * dat$maf * (1 - dat$maf) * (dat$SE^2))

print("Filtering out INDELs")
dat <- dat[!dat$INDEL, ]
print(paste("Data filtered by INDELs. Remaining Rows:", nrow(dat)))

print(paste("Filtering out INDELs. Rows:", nrow(dat), ", Columns:", ncol(dat)))

print("Filtering data based on sample size (N_cases)")
dat <- dat[dat$N_cases >= 0.7 * max(dat$N_cases), ]
print(paste("Data filtered by sample size. Remaining Rows:", nrow(dat)))

print("Filtering data based on allele frequency (EAFREQ)")
dat2 <- dat[dat$EAFREQ >= 0.005 & dat$EAFREQ <= 0.995, ]
print(paste("Data filtered by EAFREQ. Remaining Rows:", nrow(dat2)))

print("Calculating Z-scores")
dat2$z <- dat2$BETA / dat2$SE

print("Sorting data by chromosome and position")
dat2 <- dat2[order(dat2$CHR, dat2$POS), ]

print("Filtering out specific region in chromosome 11")
dat2 <- dat2[-(which(dat2$CHR == 11 & dat2$POS >= 29978453 & dat2$POS <= 80288956)), ]
print(paste("Data filtered for specific region. Remaining Rows:", nrow(dat2)))

# Further processing and tidying data
print("Tidying and processing data for output")
processed_data <- dat2 %>%
  separate(CHRBP_B37, into = c("CHR1", "BP1"), sep = ":") %>%
  transmute(
    rsID,
    CHR = CHR1, # Extract CHR from CHRBP_B37
    BP = BP1, # Extract BP from CHRBP_B37
    SNP = paste0(CHR, ":", BP1),
    A1 = EA,
    A2 = NEA,
    EAFREQ = EAFREQ,
    BETA = BETA, # No renaming needed, but included for clarity
    SE,
    N_samples = N,
    N_cases = N_cases,
    N_controls = N_controls,
    Neff = 4/((1/N_cases)+(1/N_controls)),
    P = P, # No renaming needed, but included for clarity
    NMISS = N) %>%
  filter(!is.na(CHR))

print(paste0("Original sumst SNP count: ", nrow(dat1)))
print(paste0("Processed sumst SNP count: ", nrow(processed_data)))

# Write output
print("Writing processed data to output file")
write.table(processed_data, file=opt$output, col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)
print(paste("Data written to output file successfully at", Sys.time()))