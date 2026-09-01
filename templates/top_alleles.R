#!/apps/eb/2020b/skylake/software/R/4.2.1-foss-2022a/bin/Rscript

#Install packages
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  plyr,
  RColorBrewer,
  vroom,
  tidyverse,
  qqman,
  gridExtra,
  ggplot2,
  dplyr)

pheno <-  read.table('${step_2_out}', header = T) %>% 
  #filter(str_detect(MarkerID, '^HLA')) %>%
  #select(MarkerID, p.value)
  select(CHROM, GENPOS, ID, BETA, SE, LOG10P, A1FREQ, ALLELE0, ALLELE1)

pheno\$disease <- '${pheno_label}'

# startswith <- function(string, prefix) {
#   substr(string, 1, nchar(prefix)) == prefix
# }

# pheno\$fourdgt <- ifelse(startswith(pheno\$MarkerID, 'HLA_A'), stringr::str_extract(pheno\$MarkerID, '^HLA_A\\\\*.{2}:.{2}'),
#                       ifelse(startswith(pheno\$MarkerID, 'HLA_B'), stringr::str_extract(pheno\$MarkerID, '^HLA_B\\\\*.{2}:.{2}'),
#                             ifelse(startswith(pheno\$MarkerID, 'HLA_C'), stringr::str_extract(pheno\$MarkerID, '^HLA_C\\\\*.{2}:.{2}'),
#                                 ifelse(startswith(pheno\$MarkerID, 'HLA_DQA1'), stringr::str_extract(pheno\$MarkerID, '^HLA_DQA1\\\\*.{2}:.{2}'),
#                                    ifelse(startswith(pheno\$MarkerID, 'HLA_DQB1'), stringr::str_extract(pheno\$MarkerID, '^HLA_DQB1\\\\*.{2}:.{2}'),
#                                           ifelse(startswith(pheno\$MarkerID, 'HLA_DRB1'), stringr::str_extract(pheno\$MarkerID, '^HLA_DRB1\\\\*.{2}:.{2}'),      
#                                              ifelse(startswith(pheno\$MarkerID, 'HLA_DPA1'), stringr::str_extract(pheno\$MarkerID, '^HLA_DPA1\\\\*.{2}:.{2}'),
#                                                  ifelse(startswith(pheno\$MarkerID, 'HLA_DPB1'), stringr::str_extract(pheno\$MarkerID, '^HLA_DPB1\\\\*.{2}:.{2}'),
#                                                                        '')))))))) 
pheno_nona <- na.omit(pheno)

write.table(pheno_nona, "clean_${pheno_label}.tsv", col.names = T, row.names = F, quote = F, sep = "\\t")


