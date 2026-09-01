#!/apps/eb/2020b/skylake/software/R/4.2.1-foss-2022a/bin/Rscript

#Install packages
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  readr,
  dplyr)

# Get list of file names ending with .tsv in the directory
file_names <- list.files(pattern = "\\\\.tsv\$")

# Read the first file to get column names
#first_file <- read.table(file_names[1], header = TRUE, sep = "\\t")
#column_names <- colnames(first_file)

# Read and combine files, skipping header row for subsequent files
combined_data <- lapply(file_names, read.table, header = TRUE, sep = "\\t") %>%
                 bind_rows()

# Write combined data to a new TSV file
write.table(combined_data, file = "combined_file.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)

