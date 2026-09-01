#!/apps/eb/2020b/skylake/software/R/4.2.1-foss-2022a/bin/Rscript

#Load library
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  ggplot2,
  dplyr
)

# Load in SAIGE output
saige_output <- read.table('${step_2_out}', header = TRUE)

#Define startswith function
startswith <- function(string, prefix) {
  substr(string, 1, nchar(prefix)) == prefix
}

#create a column to use for labelling based on GRCh37
saige_output\$col_Label <- ifelse(saige_output\$GENPOS >= 29942532  & saige_output\$GENPOS <= 29945870, "HLA-A",
                               ifelse(saige_output\$GENPOS >= 31353875 & saige_output\$GENPOS <= 31357179, "HLA-B",
                                      ifelse(saige_output\$GENPOS >= 31268749 & saige_output\$GENPOS <= 31272092, "HLA-C",
                                             ifelse(saige_output\$GENPOS >= 32637406 & saige_output\$GENPOS <= 32655272, "HLA-DQA1",
                                                    ifelse(saige_output\$GENPOS >= 32659467 & saige_output\$GENPOS <= 32666657, "HLA-DQB1",
                                                           ifelse(saige_output\$GENPOS >= 33064569 & saige_output\$GENPOS <= 33080748, "HLA-DPA1",
                                                                  ifelse(saige_output\$GENPOS >= 33075990 & saige_output\$GENPOS <= 33089696, "HLA-DPB1",
                                                                         ifelse(saige_output\$GENPOS >= 32578775 & saige_output\$GENPOS <= 32589848, "HLA-DRB1", "others"))))))))

#create a column to use for labelling
saige_output\$shape_Label <-ifelse(startswith(saige_output\$ID, "HLA"), "HLA",
                                ifelse(startswith(saige_output\$ID, "SNP"), "SNP",
                                       ifelse(startswith(saige_output\$ID, "AA"), "AA", "rsid")))

#Manhattan plot  

out<- ggplot(saige_output) +
  geom_point(aes(x=GENPOS,y=LOG10P, colour=col_Label, shape=shape_Label)) +
  theme_classic() +
  #geom_hline(yintercept =7.5,col="red") +
  #geom_hline(yintercept =5,col="blue") +
  scale_color_manual(values = c("#FC4E07", "#008000", "black" ,"#386CB0", "#E7B800", "blue","grey", "#00AFBB","#808080", "#CC00FFFF"))+
  ggtitle("${dataset}_${pheno_label}_regenie_Manhattan Plot")+
  labs(x = "MHC region - Chromosome 6")

pdf("${dataset}_${pheno_label}_phewas.pdf", width = 6, height = 4)
print(out)
dev.off()