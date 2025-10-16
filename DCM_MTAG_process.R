# Load libraries necessary for the script
# Adding comments to track each step of the execution in Slurm
suppressPackageStartupMessages(library(optparse))  # Parsing command-line options
suppressPackageStartupMessages(library(data.table))  # Efficient data reading and manipulation
suppressPackageStartupMessages(library(dplyr))  # Data manipulation and transformation
suppressPackageStartupMessages(library(readr))  # Reading and writing data
suppressPackageStartupMessages(library(tidyr))  # Data tidying

print("Libraries loaded successfully")

# Set up command-line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default="data/raw/DCM_MTAG_Ecc_global_LVESVi_processed.tsv", 
              help="Path to the input sumstats file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="data/sumst_processed/DCM_MTAG_sumst_37_exclMYBPC3reg.txt",
              help="Path to the output file", metavar="file")
)

print("Command-line options set up")

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))
print(paste("Input file:", opt$file))
print(paste("Output file:", opt$output))

# Set options used in the script
options(scipen = 999)
# Read data
print("Reading data from input file")
dat1 <- fread(opt$file)
print(paste("Data read successfully. Rows:", nrow(dat1), ", Columns:", ncol(dat1)))

# Process the data: select and rename columns
print("Processing data: selecting and renaming columns")
dat <- dat1

# Calculating minor allele frequency (maf)
dat$maf <- dat$EAFREQ
print("Calculating minor allele frequency (maf)")
dat[dat$EAFREQ > 0.5, 'maf'] <- 1 - dat[dat$EAFREQ > 0.5, 'maf']

print(paste("EAFREQ > 0.5. Rows:", nrow(dat), ", Columns:", ncol(dat)))

# Calculating Neff
print("Calculating Neff")
dat$Neff <- 1/(2 * dat$maf * (1 - dat$maf) * (dat$SE^2))

# Filtering data based on Neff
print("Filtering data based on Neff")
dat <- dat[dat$mtag_Neff >= 0.7 * max(dat$mtag_Neff), ]
print(paste("Data filtered by Neff. Remaining Rows:", nrow(dat)))

# Additional filtering based on EAFREQ
print("Filtering data based on EAFREQ")
dat2 <- dat[dat$EAFREQ >= 0.005 & dat$EAFREQ <= 0.995, ]
print(paste("dat$EAFREQ >= 0.005 & dat$EAFREQ <= 0.995. Rows:", nrow(dat), ", Columns:", ncol(dat)))

# Calculating Z-scores
print("Calculating Z-scores")
dat2$z <- dat2$BETA / dat2$SE

# Sorting data by chromosome and position
print("Sorting data by chromosome and position")
dat2 <- dat2[order(dat2$CHR, dat2$POS), ]

# Filtering out specific region (chromosome 11)
print("Filtering out specific region in chromosome 11")
dat2 <- dat2[-(which(dat2$CHR == 11 & dat2$POS >= 29978453 & dat2$POS <= 80288956)), ]
print(paste("Data filtered for specific region. Remaining Rows:", nrow(dat2)))

# Further processing and tidying data
print("Tidying and processing data for output")
processed_data <- dat2 %>%
  separate(CHRBP_B37, into = c("CHR1", "BP1"), sep = ":") %>%
  transmute(
    rsID,
    CHR = CHR1,
    BP = BP1,
    SNP = paste0(CHR1, ":", BP1),
    EAFREQ = EAFREQ,
    A1 = EA,
    A2 = NEA,
    BETA,
    SE,
    N_samples = NA,
    N_cases = NA,
    N_controls = NA,
    Neff = mtag_Neff,
    P) %>%
  filter(!is.na(CHR))

print(paste0("Original sumst SNP count: ", nrow(dat1)))
print(paste0("Processed sumst SNP count: ", nrow(processed_data)))

# Write output
print("Writing processed data to output file")
write.table(processed_data, file=opt$output, col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)
print("Data written to output file successfully")
print(paste("Data written to output file successfully at", Sys.time()))