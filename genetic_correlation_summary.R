######=======================================
# Load required packages
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# Define root directory with LDSC subfolders
root_dir <- "results_rg/"
subdirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)

# Function to extract genetic correlation summary line
extract_summary_block <- function(file) {
  lines <- readLines(file)
  start <- grep("Summary of Genetic Correlation Results", lines)
  end   <- grep("Analysis finished", lines)
  
  if (length(start) == 0 || length(end) == 0) return(NULL)
  
  block <- lines[(start + 1):(end[1] - 1)]
  block <- block[trimws(block) != ""]
  if (length(block) < 2) return(NULL)
  
  return(block[2])
}

# Container for master summary
all_results <- list()
# Loop through subfolders
for (subdir in subdirs) {
  log_files <- list.files(subdir, pattern = "\\.log$", full.names = TRUE)
  if (length(log_files) == 0) next
  
  summary_lines <- lapply(log_files, extract_summary_block)
  summary_lines <- summary_lines[!sapply(summary_lines, is.null)]
  if (length(summary_lines) == 0) next
  
  temp_file <- tempfile()
  writeLines(unlist(summary_lines), temp_file)
  
  df <- read.table(temp_file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                    "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                    "gcov_int", "gcov_int_se")
  
  df <- df %>%
    mutate(
      trait1 = str_extract(p1, "[^/]+(?=\\.sumstats\\.gz$)"),
      trait2 = str_extract(p2, "(?<=GWAS_sumstats_EUR__invnorm_)[^_/]+(?=__)"),

    # If trait2 is still NA, fall back to a second pattern (e.g., match last meaningful part before .tsv.gz)
    trait2 = ifelse(
      is.na(trait2),
      str_extract(p2, "(?<=/GWAS_sumstats_EUR__invnorm_)[^_/]+(?=(_TOTALsample)?\\.tsv\\.gz)"),
      trait2
    ),

    # Final fallback: try basename if still NA (for bmi_clean or other directly named files)
    trait2 = ifelse(
      is.na(trait2),
      str_extract(basename(p2), "^[^\\.]+"),
      trait2
    )
    )
  write.csv(df, paste0(subdir, "_summary.csv"), row.names = FALSE)
  cat(paste0("✅ Saved summary to: ", paste0(subdir, "_summary.csv"), "\n"))
  all_results[[basename(subdir)]] <- df
}

# List all summary CSVs from each subdirectory
#summary_files <- list.files("results_rg", pattern = "_summary.csv$", full.names = TRUE, recursive = TRUE)
#
# Initialize combined data frame
#all_summaries <- purrr::map_dfr(summary_files, read_csv)
all_summaries <- imap_dfr(all_results, ~ {
  df <- .x
  df$p1 <- .y  # Add name of list element as trait1
  df
})

all_summaries %>% filter(is.na(trait2))

# Summarise duplicates (mean if repeated runs)
all_summaries <- all_summaries %>%
  group_by(trait2, trait1) %>%
  summarise(rg = mean(rg, na.rm = TRUE),
            p = mean(p, na.rm = TRUE),
            .groups = "drop")

# Pivot to wide format
wide_summary <- all_summaries %>%
  pivot_wider(
    names_from = trait1,
    values_from = c(rg, p),
    names_glue = "{.value}_{trait1}",
    values_fill = NA  # 🔑 Ensures missing pairs are retained as NA
  )

library(readr)
# Write overall summary
write_csv(wide_summary, "eur_ref_genetic_correlation_traits.csv")

cat("✅ Saved combined summary table to 'genetic_correlation_traits.csv'\n")


######======================================= results_rg_pairwise
# Load required packages
library(dplyr)
library(stringr)

# Define root directory with LDSC subfolders
root_dir <- "results_rg_pairwise/"

# Function to extract genetic correlation summary line
extract_summary_block <- function(file) {
  lines <- readLines(file)
  start <- grep("Summary of Genetic Correlation Results", lines)
  end   <- grep("Analysis finished", lines)
  
  if (length(start) == 0 || length(end) == 0) return(NULL)
  
  block <- lines[(start + 1):(end[1] - 1)]
  block <- block[trimws(block) != ""]
  if (length(block) < 2) return(NULL)
  
  return(block[2])
}

# Container for master summary
all_results <- list()
# Loop through subfolders
log_files <- list.files(root_dir, pattern = "^eur_ref.*\\.log$", full.names = TRUE)
# log_files <- list.files(root_dir, pattern = "\\.log$", full.names = TRUE)

if (length(log_files) == 0) next

summary_lines <- lapply(log_files, extract_summary_block)
summary_lines <- summary_lines[!sapply(summary_lines, is.null)]
if (length(summary_lines) == 0) next

temp_file <- tempfile()
writeLines(unlist(summary_lines), temp_file)

df <- read.table(temp_file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                  "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                  "gcov_int", "gcov_int_se")

df <- df %>%
  mutate(
    trait1 = str_extract(p1, "[^/]+(?=\\.sumstats\\.gz$)"),
    trait2 = str_extract(p2, "[^/]+(?=\\.sumstats\\.gz$)"),
  )
write.csv(df, paste0(root_dir,"eur_ref_results_rg_pairwise_summary.csv"), row.names = FALSE)
cat(paste0("✅ Saved summary to: ", paste0(root_dir,"results_rg_pairwise_summary.csv"), "\n"))

######======================================= results_rg_gwas
# Load required packages
library(dplyr)
library(stringr)

# Define root directory with LDSC subfolders
root_dir <- "results_rg_gwas/"

# Function to extract genetic correlation summary line
extract_summary_block <- function(file) {
  lines <- readLines(file)
  start <- grep("Summary of Genetic Correlation Results", lines)
  end   <- grep("Analysis finished", lines)
  
  if (length(start) == 0 || length(end) == 0) return(NULL)
  
  block <- lines[(start + 1):(end[1] - 1)]
  block <- block[trimws(block) != ""]
  if (length(block) < 2) return(NULL)
  
  return(block[2])
}

# Container for master summary
all_results <- list()
# Loop through subfolders
log_files <- list.files(root_dir, full.names = TRUE)
# log_files <- list.files(root_dir, pattern = "\\.log$", full.names = TRUE)

if (length(log_files) == 0) next

summary_lines <- lapply(log_files, extract_summary_block)
summary_lines <- summary_lines[!sapply(summary_lines, is.null)]
if (length(summary_lines) == 0) next

temp_file <- tempfile()
writeLines(unlist(summary_lines), temp_file)

df <- read.table(temp_file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                  "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                  "gcov_int", "gcov_int_se")

df <- df %>%
  mutate(
    trait1 = str_extract(p1, "[^/]+(?=\\.sumstats\\.gz$)"),
    trait2 = str_extract(p2, "[^/]+(?=\\.sumstats\\.gz$)"),
  )
write.csv(df, paste0(root_dir,"LV_traits_pairwise_summary.csv"), row.names = FALSE)
cat(paste0("✅ Saved summary to: ", paste0(root_dir,"LV_traits_pairwise_summary.csv"), "\n"))

