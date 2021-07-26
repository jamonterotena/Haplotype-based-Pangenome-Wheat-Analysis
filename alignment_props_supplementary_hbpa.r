## RUN THE SCRIPT FUNCTIONS

#################################################################################################################################################### 
############################################################### DATA GENERATION ####################################################################
#################################################################################################################################################### 
chromosome <- c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D")
chromosome_column <- character()
for(j in 1:length(chromosome)){
  chromosome_column <- append(chromosome_column, rep(paste0("chr", chromosome[[j]], collapse = ""), 20))
}

comparison <- c("CLA-CLA", "CLA-SLA")
comparison_column <- character()
for(k in 1:length(comparison)){
  comparison_column <- append(comparison_column, rep(comparison[[k]], 10))
}
comparison_column <- rep(comparison_column, 21)

cultivar <- c("arinalrfor", "jagger", "julius", "lancer", "landmark", "mace", "norin61", "stanley", "sy_mattis", "chinese")
cultivar_column <- rep(cultivar, 42)

length(comparison_column) == length(chromosome_column)
length(chromosome_column) == length(cultivar_column)

perc_cov_column <- numeric()
length_aln_column <- numeric()
number_aln_column <- numeric()

table_cultivar_chr_length <- read.table(file = "table_cultivar_chr_length.txt", header = TRUE)
path <- c("C:/Users/Windows 10/Documents/RStudio-thesis/all_chr_rds/all_20_kb_filtered_delta_CHROMOSOME_tables.rds")

chromosome <- c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D")

for (i in 1:length(chromosome)){
  print(paste0("chr", chromosome[i], ". Reading the RDS file:", sep = ""))
  data <- readRDS(file = gsub("CHROMOSOME", chromosome[i], path))
  print("Done")
  print(paste0("chr", chromosome[i], ". Calculating average number of alignments:", sep = ""))
  output.3 <- calculate_av_number_aln(data = data)
  print("Done")
  print(paste0("chr", chromosome[i], ". Calculating average percentage of coverage by alignments:", sep = ""))
  output.1 <- calculate_perc_coverage(data = data, chrlength_v_cultivar = table_cultivar_chr_length, target_chr = paste0("chr", chromosome[i], collapse = ""))
  print("Done")
  print(paste0("chr", chromosome[i], ". Calculating average alignment length:", sep = ""))
  output.2 <- calculate_av_aln_length(data = data)
  print("Done")
  perc_cov_column <- append(perc_cov_column, as.numeric(output.1[[2]]))
  perc_cov_column <- append(perc_cov_column, as.numeric(output.1[[4]]))
  length_aln_column <- append(length_aln_column, as.numeric(output.2[[2]]))
  length_aln_column <- append(length_aln_column, as.numeric(output.2[[4]]))
  number_aln_column <- append(number_aln_column, as.numeric(output.3[[2]]))
  number_aln_column <- append(number_aln_column, as.numeric(output.3[[4]]))
}

perc_cov_aln_table <- data.frame("chromosome" = chromosome_column, "comparison_type" = comparison_column, "reference_cultivar_name" = cultivar_column, "perc_coverage_ref_genome_aln" = perc_cov_column)
length_aln_table <- data.frame("chromosome" = chromosome_column, "comparison_type" = comparison_column, "reference_cultivar_name" = cultivar_column, "av_aln_length" = length_aln_column)
number_aln_table <- data.frame("chromosome" = chromosome_column, "comparison_type" = comparison_column, "reference_cultivar_name" = cultivar_column, "av_aln_number" = number_aln_column)
aln_prop_full_table <- data.frame("chromosome" = chromosome_column, "comparison_type" = comparison_column, "reference_cultivar_name" = cultivar_column, "perc_coverage_ref_genome_aln" = perc_cov_column, "av_aln_length" = length_aln_column, "av_aln_number" = number_aln_column)

write.csv2(x = perc_cov_aln_table, file = "perc_cov_aln_table.csv", row.names = FALSE)
write.csv2(x = length_aln_table, file = "length_aln_table.csv", row.names = FALSE)
write.csv2(x = number_aln_table, file = "number_aln_table.csv", row.names = FALSE)
write.csv2(x = aln_prop_full_table, file = "aln_prop_full_table.csv", row.names = FALSE)

#################################################################################################################################################### 
############################################################### STATISTICAL ANALYSIS ###############################################################
#################################################################################################################################################### 

library(dplyr)
library(magrittr)
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(viridis)
library(stringr)
library(multcomp)
library(ggpubr)

############################################################### DATA PREP ############################################################### 
perc_cov_raw <- read.table(file = "perc_cov_aln_table.csv", header = TRUE, sep = ";")
str(perc_cov_raw)
values <- gsub(",", ".", perc_cov_raw$perc_coverage_ref_genome_aln)
values <- as.numeric(values)
perc_cov_raw$perc_coverage_ref_genome_aln <- NA
perc_cov_raw$perc_coverage_ref_genome_aln <- values

############################################################### BOXPLOTS ############################################################### 
par(mfrow=c(1,1))
# CLA-CLA vs CLA-SLA (N=210, data from references)
boxplot ( perc_coverage_ref_genome_aln ~ comparison_type, range = 0, data = perc_cov_raw)
# CLA-CLA vs CLA-SLA (N=21, data from chromosomes)
CHR <- unique(perc_cov_raw$chromosome)
by_chr <- numeric()
for (i in 1:length(CHR)){
  by_chr <- append(by_chr, mean(subset(perc_cov_raw$perc_coverage_ref_genome_aln, perc_cov_raw$comparison_type == "CLA-CLA" & perc_cov_raw$chromosome == CHR[i])))
  by_chr <- append(by_chr, mean(subset(perc_cov_raw$perc_coverage_ref_genome_aln, perc_cov_raw$comparison_type == "CLA-SLA" & perc_cov_raw$chromosome == CHR[i])))
                   
}
labels <- rep(c("vs CHROMOSOME-LEVEL ASSEMBLIES", "vs SCAFFOLD-LEVEL ASSEMBLIES"), length(CHR))
perc_cov_by_chr <- data.frame (comparison = labels, coverage = by_chr)
g <- ggplot(perc_cov_by_chr, aes(x=as.factor(comparison), y=coverage)) + 
  geom_boxplot(color = "black", fill = "lightgrey", alpha=0.8, outlier.shape = TRUE) +
  ylim(0,80) +
  ylab("Percentage of alignment coverage of the reference") +
  xlab("Query assembly in pairwise alignments") +
  ggtitle(label = "AVERAGE PERCENTAGE OF ALIGNMENT COVERAGE OF THE REFERENCE CHROMOSOME") +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 0.02) +
  theme(plot.title   = element_text(size = 18, face = "bold.italic"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"))
#panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "black"), 
#panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "black"))
g
my_comparisons <- list(c("vs CHROMOSOME-LEVEL ASSEMBLIES", "vs SCAFFOLD-LEVEL ASSEMBLIES"))
g <- g + stat_compare_means(comparisons = my_comparisons, label.y = 75, label = "p.signif")
g
#CLA-CLA chr5B vs all
perc_cov_only_CLA_CLA <- perc_cov_raw[!grepl("CLA-SLA", perc_cov_raw$comparison_type),]
bymedian <- with(perc_cov_only_CLA_CLA, reorder(perc_cov_only_CLA_CLA$chromosome, perc_cov_only_CLA_CLA$perc_coverage_ref_genome_aln, median, na.rm = TRUE))
colnames(perc_cov_only_CLA_CLA) <- c("chromosome", "comparison", "reference", "coverage")
p <- ggplot(perc_cov_only_CLA_CLA, aes(x=as.factor(bymedian), y=coverage)) + 
  geom_boxplot(color = "black", fill = "lightgrey", alpha=0.8, outlier.shape = TRUE) +
  ylim(50,80) +
  ylab("Percentage of alignment coverage of the reference") +
  xlab("Reference chromosome") +
  ggtitle(label = "AVERAGE PERCENTAGE OF ALIGNMENT COVERAGE OF THE REFERENCE CHROMOSOME IN COMPARISONS BETWEEN CHR-LEVEL ASSEMBLIES") +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 0.5) +
  theme(plot.title   = element_text(size = 15, face = "bold.italic"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'dashed', colour = "grey"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'dashed', colour = "grey"))
my_comparisons <- list(c("vs CHROMOSOME-LEVEL ASSEMBLIES", "vs SCAFFOLD-LEVEL ASSEMBLIES"))
p <- p + annotate ("text", x = 1, y = 56, label = "*")
p
#CLA-SLA
#perc_cov_only_CLA_SLA <- perc_cov_raw[!grepl("CLA-CLA", perc_cov_raw$comparison_type),]
#bymedian <- with(perc_cov_only_CLA_SLA, reorder(perc_cov_only_CLA_CLA$chromosome, perc_cov_only_CLA_SLA$perc_coverage_ref_genome_aln, median, na.rm = TRUE))
#boxplot ( perc_coverage_ref_genome_aln ~ bymedian, data = perc_cov_only_CLA_SLA, range = 0)
#
#aspect_ratio <- 2.5
#height <- 7
#ggsave(g, height = 7 , width = 7 * aspect_ratio)
############################################################### COMPARING BETWEEN COMPARISONS ###############################################################
# Descriptive statistics
str(perc_cov_by_chr)
perc_cov_by_chr <- perc_cov_by_chr[complete.cases(perc_cov_by_chr), ]
perc_cov_by_chr %>% group_by(comparison) %>%
  summarise (n=n(), mean=mean(coverage), var=var(coverage),
             min=min(coverage), max=max(coverage))

attach(perc_cov_by_chr)
## Checking for prerequisites
# FIRST WAY: SINGLE-FACTOR ANOVA (CHECKING FOR PREREQUISITES)
bartlett.test(coverage ~ comparison) # The variance of the mean is not homogeneous (heteroscedasticity)
shapiro.test(perc_cov_by_chr$coverage[grepl ("vs CHROMOSOME-LEVEL ASSEMBLIES", perc_cov_by_chr$comparison)])
shapiro.test(perc_cov_by_chr$coverage[grepl ("vs SCAFFOLD-LEVEL ASSEMBLIES", perc_cov_by_chr$comparison)]) # The sample is not normally distributed
wilcox.test(coverage ~ comparison) # Heteroscedasticity and not-normal distribution. Mean coverage is significantly different (p-value 0.05) between comparisons

# SECOND WAY: ANALYSIS OF RESIDUALS FOR AN ANOVA (CHECKING FOR PREREQUISITES)
resid <- residuals ( aov (coverage ~ comparison) )
shapiro.test ( resid ) # Data not normally distributed
bartlett.test ( resid ~ comparison ) # Heteroscedastic data
par(mfrow=c(2,2))
plot( aov ( coverage ~ comparison ) )
# Prerequisites for ANOVA test are not fulfilled, so we are going to use the wilcox test to study differences between the two samples

# Statistical test
wilcox.test(coverage ~ comparison) # Heteroscedasticity and not-normal distribution. Mean coverage is significantly different (p-value 0.05) between comparisons
# kruskal.test(coverage ~ comparison) # No assumption on distribution or variances. Mean coverage is significantly different (p-value 0.05) between comparisons
detach(perc_cov_by_chr)

############################################################### COMPARING BETWEEN CHROMOSOMES ###############################################################
# Descriptive statistics
colnames(perc_cov_only_CLA_CLA) <- c("chromosome", "comparison", "reference", "coverage")
str(perc_cov_only_CLA_CLA)
perc_cov_only_CLA_CLA <- perc_cov_only_CLA_CLA[complete.cases(perc_cov_only_CLA_CLA), ]
sum <- perc_cov_only_CLA_CLA %>% group_by(chromosome) %>%
  summarise (n=n(), mean=mean(coverage), var=var(coverage),
             min=min(coverage), max=max(coverage))
sum # Mean values, variances, mins and maxs of coverage in each chromosome
attach(perc_cov_only_CLA_CLA)

## Checking for prerequisites
# FIRST WAY: SINGLE-FACTOR ANOVA (CHECKING FOR PREREQUISITES)

# SECOND WAY: ANALYSIS OF RESIDUALS FOR AN ANOVA (CHECKING FOR PREREQUISITES)
resid <- residuals ( aov (coverage ~ chromosome) )
shapiro.test ( resid ) # Data not normally distributed
bartlett.test ( resid ~ chromosome ) # Heteroscedastic data
par(mfrow=c(2,2))
plot( aov ( coverage ~ chromosome ) )
# Prerequisites for ANOVA test are not fulfilled, so we are going to use the non-parametric Kruskal test and the pairwise wilcox test to study differences in cov between chromosomes

# Statistical test
kruskal.test(coverage ~ chromosome) # Mean coverage is significantly different (p-value 0.05) in at least one of the chromosomes (no assumption on distribution or variances)
pairwise.wilcox.test( coverage , chromosome ) # Shows differences between chromosomes. Chr5B appears to be different to the rest of the chromosomes (p-value 0.0325)
############## 
detach(perc_cov_only_CLA_CLA)
