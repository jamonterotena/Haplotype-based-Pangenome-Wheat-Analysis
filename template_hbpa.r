########################################        HAPLOTYPE-BASED PANGENOME ANALYSIS IN WHEAT (template)        ########################################         


# As the writer of this script I would like to express appreciation and give all credit to the authors of the following article for creating original raw data and scripts, which are required for this R script.
# Brinton, J., Ramirez-Gonzalez, R. H., Simmonds, J., Wingen, L., Orford, S., Griffiths, S., 10 Wheat Genome Project, Haberer, G., Spannagl, M., Walkowiak, S., Pozniak, C., & Uauy, C. (2020). A haplotype-led approach to increase the precision of wheat breeding. Communications biology, 3(1), 712. https://doi.org/10.1038/s42003-020-01413-2
# This script was written in R 4.1.0 by Jose Antonio Montero Tena as part of his master thesis in 2021. Free use and modification of this script for personal interest is encouraged.



# INTRODUCTION: SNP-haplotypes discovered from marker arrays can be redefined as IBS-haplotypes, which are identical-by-state sequences sharing 100 percentage of identity between the genomes and discovered by sequence comparison. Brinton et al developed a method in 2020 that defined haploblocks as physical regions with >= 99.99 pident across fixed-size bins of 5-, 2.5- or 1-Mbp between pairwise comparisons of cultivars in the 10+ Wheat Genome Project. This method will subsequently be referred to as the Brinton's method. In these IBS-haploblocks, IBS-haplotypes can be either shared haplotypes, when pairs of assemblies share one identical-by-state sequence in this region, or unique haplotypes, when no other cultivar presents the same sequence in the haploblock. Additionally, Brinton and colleagues published crop-haplotypes.com, a website that provides an interactive graphic visualization of the shared haplotypes between the wheat genome assemblies. These assemblies were chromosome-level genome assemblies of nine wheat cultivars (ArinaLrFor, Jagger, Julius, Lancer, Landmark, Mace, Norin61, Stanley, SY-Mattis) and the Chinese Spring RefSev1.0 assembly alongside scaffold-level assemblies of five additional cultivars (Cadenza, Claire, Paragon, Robigus and Weebill).
# AIM: The aim of this script is to use Brinton's method and resources to conduct both a chromosome-scale and small-scale analysis of the pairwise comparisons between cultivars. The small-scale analysis is an original method to this script and allows the analysis of target regions within chromosomes. The aim of the small-scale analysis is to map physical start and end coordinates of the predicted haploblock both in the reference and query genomes by considering the region coverage with alignments, in other words, what is the number of alignments the haplotype prediction is based on and if there is any discontinuity that could cause false positive predictions. Also, this script aims to identify the genes contained in mapped haploblock, both in the haplotype-carrying genomes and in the IWGSC Chinese Spring Ref v1.1. Brinton et al called haplotypes both from mummer and gene-based BLAST pairwise alignments and selected for matching positions the longer predictions. Both types of analysis are provided in this script, although the results from these methods do not always match.

##### RUNNING THE SCRIPT #####  

# The script is designed so that you only have to edit only the upcoming parameters. You can run the script by lines (recommended for first time) or simply pressing 'Source' in the upper right corner of this window. The running time is expected to be between 5 and 7 minutes. Please, be patient and only look for solutions if error messages interrupt the pipeline. If errors do not stop coming up, you can contact the script editor at the e-mail adress jmonterotena@gmail.com. 

##### DEFINING PARAMETERS ##### 

chromosome <- "5B" # Copy the same format
reference_cultivar <- "Lancer" # Alignment coordinates will apply for this cultivar's chromosome. Can NEVER be a scaffold-level assembly (Cadenza, Claire, Paragon, Robigus, Weebil). Follow the next naming guidelines: "arinalrfor", "jagger", "julius", "lancer", "landmark", "mace", "norin61", "stanley", "sy_mattis", "chinese" (works with capital letters)
query_cultivar <- "Paragon" # This genome is aligned against the reference chromosome. Can be a scaffold-level assembly. Follow the next naming guidelines: "cadenza", "claire", "paragon", "robigus", "weebil" (works with capital letters)
zoom_start <-  655000000 # Start of the region highlighted in the small-scale analysis, in which haploblocks are sought for. ADVICE: use multiples of 5 so that the bins match with those on crop-haplotypes
zoom_end <- 665000000 # End of the region highlighted in the small-scale analysis, in which haploblocks are sought for. ADVICE: use multiples of 5 so that the bins match with those on crop-haplotypes
target_start <- 655700000 # If you have evidence of a haplotype within your zoom region that you want to compare with the IBS-haplotypes obtained in this script, write its start coordinate in the reference chromosome. If not, simply write NA. The target region must be within the zoom region
target_end <- 656600000 # If you have evidence of a haplotype within your zoom region that you want to compare with the IBS-haplotypes obtained in this script, write its end coordinate in the reference chromosome. If not, simply write NA. The target region must be within the zoom region

##### PREPARING FILES #####  

# - Script 'functions.hbpa.r' (downloadable from https://github.com/MonteroJLU/Haplotype-based-Pangenome-Wheat-Analysis.git) 
# - Raw delta files in 'all_20_kb_filtered_delta_CHROMOSOME_tables.rds'; in this case the chromosome is 5B (downloadable from https://opendata.earlham.ac.uk/wheat/under_license/toronto/Brinton_etal_2020-05-20-Haplotypes-for-wheat-breeding/nucmer/). Brinton et al filtered out alignments under 20000Kbp of length in order to get rid of non-syntenic retrotransposons
# - Text file 'table_cultivar_chr_length.txt' (downloadable from https://github.com/MonteroJLU/Haplotype-based-Pangenome-Wheat-Analysis.git)
# - File 'projectedGenes__Triticum_aestivum_REFERENCECULTIVAR_v1.0.gff', in this case the reference is LongReach Lancer (downloadable from https://webblast.ipk-gatersleben.de/downloads/wheat/gene_projection/)
# - File 'geneid_2_chinese.sourceid.txt' (downloadable from https://webblast.ipk-gatersleben.de/downloads/wheat/gene_projection/)
# - File 'varieties_all_identities_2000bp.tar.gz' (downloadable from https://opendata.earlham.ac.uk/wheat/under_license/toronto/Brinton_etal_2020-05-20-Haplotypes-for-wheat-breeding/pairwise_blast/)
# - File 'iwgsc_refseq_v1.2_gene_annotation.zip' (downloadable from https://urgi.versailles.inrae.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.2/)

##### RUNNING FUNCTIONS #####  

# Source the script 'functions.hbpa.r' to read the packages and functions required for this script.

############### 1. CALLING HAPLOTYPES FROM RAW DATA CONTAINING >= 20-KBP-LONG MUMMER PAIRWISE ALIGNMENTS ###############   

print("1. CALLING HAPLOTYPES FROM RAW DATA CONTAINING >= 20-KBP-LONG MUMMER PAIRWISE ALIGNMENTS")

##########  1.1. CHROMOSOME-SCALE ANALYSIS ##########  

print("1.1. CHROMOSOME-SCALE ANALYSIS")
print("**1.1.1. Download and save the required documents in the working directory**")
print("1.1.2. Run the script 'functions_hbpa.r'")
print("1.1.3. Read the rds file into a data frame")
aln_library <- readRDS(file = paste0("all_20_kb_filtered_delta_", chromosome,"_tables.rds"))
print("1.1.4. Create a subset of the data frame containing only the alignments from our pairwise comparison")
aln_subset <- aln_library[grepl(paste0("^", tolower(reference_cultivar), sep = ""), aln_library$comparison) & grepl(tolower(query_cultivar), aln_library$comparison),]
print("1.1.5. Scatter-plot the alignment midpoints (X: r_mid, Y: q_mid). Check plots")
plot_diagonal_scatterplot(aln_subset, cap_lower = 90.00, cap_upper = 100, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 100000000)
print("1.1.6. Dot-plot the alignments to show percentage of identity and alignment length (X: r_mid, Y: perc_id). Check plots")
plot_aln_pid_and_length(aln_subset, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 100000000, dot_size = 2)
print("1.1.7. Check for alignment properties in the data frame (percentage of alignment coverage of reference chromosome, average aln length and av aln number)")
table_chrlength_v_cultivar <- read.table(file = "table_cultivar_chr_length.txt", sep = "\t", header = TRUE)
chr_length <- table_chrlength_v_cultivar$sequence_length[grepl(tolower(reference_cultivar), table_chrlength_v_cultivar$cultivar_name) & grepl(paste0("chr", unique(aln_subset$chrom)), table_chrlength_v_cultivar$sequence_name)]
print(paste0("chr", unique(aln_subset$chrom), " is ", chr_length, " bp-long in ", reference_cultivar))
print(paste0(round(mean(aln_subset$r_length), 0), " is the average alignment length in ", reference_cultivar, "-", query_cultivar, paste0(" chr", unique(aln_subset$chrom)), " comparison"))
print(paste0(nrow(aln_subset), " is the number of alignments in ", reference_cultivar, "-", query_cultivar, paste0(" chr", unique(aln_subset$chrom)), " comparison")) 
print(paste0(round((sum(aln_subset$r_length)/chr_length*100), 0), "% is the alignment coverage in ", reference_cultivar, "-", query_cultivar, paste0(" chr", unique(aln_subset$chrom)), " comparison")) 
coverage <- COV(aln_library, chr_length = table_chrlength_v_cultivar)
print(coverage)
length <- LENGTH(aln_library)
print(length)
number <- NUMBER(aln_library)
print(number)
bin_size <- c(5000000, 2500000, 1000000)
names(bin_size) <- c("bin size: 5-Mbp", "bin size: 2.5-Mbp", "bin size: 1-Mbp")
for (i in 1:3){
  print("average expected number of alignments per bin across the chromosome")
  print(round(nrow(aln_subset)/(chr_length/bin_size[i]), 1))
}
print("1.1.8. Boxplot the alignments (X: bin, Y: perc_id_median). Check plots")
plot_boxplots_bin_median(aln_subset, bin_size = 10000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 50000000, show_outliers = FALSE)
plot_boxplots_bin_median(aln_subset, bin_size = 5000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 50000000, show_outliers = FALSE)
plot_boxplots_bin_median(aln_subset, bin_size = 2500000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 50000000, show_outliers = FALSE)
plot_boxplots_bin_median(aln_subset, bin_size = 1000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 50000000, show_outliers = FALSE)
print("1.1.9. Plot the median percentage of identity across the chromosome or median line (X: r_mid, Y: perc_id_median)")
medians_aln_subset_10Mbp <- plot_line_bin_median(aln_subset, bin_size = 10000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 50000000)
medians_aln_subset_10Mbp
medians_aln_subset_5Mbp <- plot_line_bin_median(aln_subset, bin_size = 5000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 50000000)
medians_aln_subset_5Mbp
medians_aln_subset_2.5Mbp <- plot_line_bin_median(aln_subset, bin_size = 2500000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 50000000)
medians_aln_subset_2.5Mbp
medians_aln_subset_1Mbp <- plot_line_bin_median(aln_subset, bin_size = 1000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 50000000)
medians_aln_subset_1Mbp
print("1.1.10. Print information about the haploblock predictions in the pairwise comparison across the reference chromosome")
medians_aln_subset_10Mbp_bin_info <- assign_blocks_mummer(median_cutoffs = medians_aln_subset_10Mbp[["data"]], original_file = aln_subset)
medians_aln_subset_5Mbp_bin_info <- assign_blocks_mummer(median_cutoffs = medians_aln_subset_5Mbp[["data"]], original_file = aln_subset)
medians_aln_subset_2.5Mbp_bin_info <- assign_blocks_mummer(median_cutoffs = medians_aln_subset_2.5Mbp[["data"]], original_file = aln_subset)
medians_aln_subset_1Mbp_bin_info <- assign_blocks_mummer(median_cutoffs = medians_aln_subset_1Mbp[["data"]], original_file = aln_subset)
medians_aln_subset_10Mbp_block_summary <- block_summary(medians_aln_subset_10Mbp_bin_info, bin_size = 10000000, reference_name = reference_cultivar, query_name = query_cultivar )
medians_aln_subset_5Mbp_block_summary <- block_summary(medians_aln_subset_5Mbp_bin_info, bin_size = 5000000, reference_name = reference_cultivar, query_name = query_cultivar )
medians_aln_subset_2.5Mbp_block_summary <- block_summary(medians_aln_subset_2.5Mbp_bin_info, bin_size = 2500000, reference_name = reference_cultivar, query_name = query_cultivar )
medians_aln_subset_1Mbp_block_summary <- block_summary(medians_aln_subset_1Mbp_bin_info, bin_size = 1000000, reference_name = reference_cultivar, query_name = query_cultivar )
print(medians_aln_subset_10Mbp_block_summary)
print(medians_aln_subset_5Mbp_block_summary)
print(medians_aln_subset_2.5Mbp_block_summary)
print(medians_aln_subset_1Mbp_block_summary)


##########  1.2. SMALL-SCALE ANALYSIS ##########  


print("1.2. SMALL-SCALE ANALYSIS")
print("1.2.1. Scatter-plot the alignment midpoints across the zoom region (X: r_mid, Y: q_mid). Check plots")
plot_diagonal_scatterplot(aln_subset, xmin = zoom_start, xmax = zoom_end, cap_lower = 90.00, cap_upper = 100, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 1000000)
print("1.2.2. Dot-plot the alignments to show percentage of identity and alignment length in the zoom region (X: r_mid, Y: perc_id). Check plots")
graph <- plot_aln_pid_and_length(data = aln_subset, xmin = zoom_start, xmax = zoom_end, ymin = 98.5, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 1000000, dot_size = 4)
graph
print("1.2.3. Check for alignment properties in the zoom region")
aln_target <- aln_subset[(aln_subset$r_mid >= target_start) & (aln_subset$r_mid <= target_end),]
print(paste0(round(mean(aln_target$r_length), 0), " is the average alignment length for the target region between ", target_start/1e06, " and ", target_end/1e06, " Mbp in ", reference_cultivar, "-", query_cultivar, paste0(" chr", unique(aln_target$chrom)), " comparison")) 
print(paste0(nrow(aln_target), " is the number of alignments for the target region between ", target_start/1e06, " and ", target_end/1e06, " Mbp in ", reference_cultivar, "-", query_cultivar, paste0(" chr", unique(aln_target$chrom)), " comparison"))
print(paste0(round((sum(aln_target$r_length)/(target_end-target_start)*100), 0), "% is the alignment coverage for the target region between ", target_start/1e06, " and ", target_end/1e06, " Mbp in ", reference_cultivar, "-", query_cultivar, paste0(" chr", unique(aln_target$chrom)), " comparison") )
aln_zoom <- aln_subset[(aln_subset$r_mid >= zoom_start) & (aln_subset$r_mid <= zoom_end),]
print(paste0(round(mean(aln_zoom$r_length), 0), " is the average alignment length for the zoom region between ", zoom_start/1e06, " and ", zoom_end/1e06, " Mbp in ", reference_cultivar, "-", query_cultivar, paste0(" chr", unique(aln_subset$chrom)), " comparison")) 
print(paste0(nrow(aln_zoom), " is the number of alignments for the zoom region between ", zoom_start/1e06, " and ", zoom_end/1e06, " Mbp in ", reference_cultivar, "-", query_cultivar, paste0(" chr", unique(aln_subset$chrom)), " comparison"))
print(paste0(round((sum(aln_zoom$r_length)/(zoom_end-zoom_start)*100), 0), "% is the alignment coverage for the zoom region between ", zoom_start/1e06, " and ", zoom_end/1e06, " Mbp in ", reference_cultivar, "-", query_cultivar, paste0(" chr", unique(aln_subset$chrom)), " comparison")) 
bin_size <- c(5000000, 2500000, 1000000)
names(bin_size) <- c("bin size: 5-Mbp", "bin size: 2.5-Mbp", "bin size: 1-Mbp")
for (i in 1:3){
  print("average expected number of alignments per bin across zoom region")
  print(nrow(aln_zoom)/((zoom_end-zoom_start)/bin_size[i]))
}
print("1.2.4. Boxplot the bin median percentage of identity in the zoom region to check for outliers. Check plots")
plot_boxplots_bin_median(aln_subset, bin_size = 5000000, bin_start = zoom_start, bin_end = zoom_end, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 1000000, show_outliers = TRUE)
plot_boxplots_bin_median(aln_subset, bin_size = 2500000, bin_start = zoom_start, bin_end = zoom_end, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 1000000, show_outliers = TRUE)
plot_boxplots_bin_median(aln_subset, bin_size = 1000000, bin_start = zoom_start, bin_end = zoom_end, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 1000000, show_outliers = TRUE)
print("1.2.5. Plot the median line in the zoom region to see the haploblock predictions at different bin sizes. Check plots")
plot_line_bin_median(aln_subset, bin_size = 5000000, bin_start = zoom_start, bin_end = zoom_end, ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 1000000)
plot_line_bin_median(aln_subset, bin_size = 2500000, bin_start = zoom_start, bin_end = zoom_end, ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 1000000)
plot_line_bin_median(aln_subset, bin_size = 1000000, bin_start = zoom_start, bin_end = zoom_end, ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar , x_label_gap = 1000000)
print("1.2.6. Dot-plot the zoom region with the haploblock predictions and the target region and make decisions regarding the start and end limits of the IBS-haploblock. Check plots")
target <- data.frame(target_start, target_end)
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 5000000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar, dot_size = 4, x_label_gap = 1000000)
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 2500000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar, dot_size = 4, x_label_gap = 1000000)
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 1000000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar, dot_size = 4, x_label_gap = 1000000)

##### DEFINING NEW PARAMETERS ##### 

selected_start <-  # Genes will be extracted from this position. If you have no interest in redefining your region, simply write 'target_start' or 'zoom_start', to keep with the previous coordinates
selected_end <-  # Genes will be extracted until this position. If you have no interest in redefining your region, simply write 'target_end' or 'zoom_end', to keep with the previous coordinates
target_text <- "" # Text to print on the selected region

selected_haploblock <- data.frame(selected_start, selected_end)
print(selected_haploblock)
print("1.2.7. Plot summary graphs. Check plots")
plot_aln_and_bins(print_tables = FALSE, aln_subset = aln_subset, bin_size = 1000000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = selected_haploblock, target_text = target_text, fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar, dot_size = 4, x_label_gap = 1000000)
plot_bins_and_selected_region(print_tables = FALSE, aln_subset = aln_subset, bin_size = 1000000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = selected_haploblock, target_text = target_text, fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar, dot_size = 4, x_label_gap = 1000000)
plot_aln_and_bins(print_tables = FALSE, aln_subset = aln_subset, bin_size = 5000000, zoom_start = 0, zoom_end = max(aln_subset$re), highlighted_target = selected_haploblock, target_text = target_text, fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar, dot_size = 4, x_label_gap = 100000000, prediction_text = FALSE, aln_text = F)
plot_bins_and_selected_region(print_tables = FALSE, aln_subset = aln_subset, bin_size = 5000000, zoom_start = 0, zoom_end = max(aln_subset$re), highlighted_target = selected_haploblock, target_text = target_text, fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_cultivar, query_name = query_cultivar, dot_size = 4, x_label_gap = 100000000, prediction_text = FALSE)


############### 2. IDENTIFY GENES WITHIN THE JUST-MAPPED HAPLOBLOCK ###############   


print("2. IDENTIFY GENES WITHIN THE JUST-MAPPED HAPLOBLOCK")   
print("2.1. Download the files and save them in the working directory")
print("2.2. Read the gff file containing the gene model projection for the reference cultivar")
ref_gff <- read.table ( file = "projectedGenes__Triticum_aestivum_LongReach_Lancer_v1.0.gff", sep = "\t" , header = F, stringsAsFactors = F)
print("2.3. Create subsets for gene projections only and edit the table")
ref_gff_gene_only <- ref_gff[ref_gff$V3 == "gene",]
colnames(ref_gff_gene_only) <- c("chr", "annotation", "biotype", "start", "end", "score", "strand", "info", "ref_id")
ref_gff_gene_only$var_id <- gsub("ID=", "", ref_gff_gene_only$ref_id)
print("2.4. Create a subset for the haploblock") 
my_genes <- ref_gff_gene_only[grep(chromosome, ref_gff_gene_only$chr),]
my_genes <- my_genes[(my_genes$start >= selected_start) & (my_genes$end <= selected_end),]
print("2.5. Download 'geneid_2_chinese.sourceid.txt' from the last link and read the table") 
cs_id <- read.table(file = "geneid_2_chinese.sourceid.txt", sep="\t", header=T, stringsAsFactors = F)
print("2.6. Add a column to the subset with the IDs of the gene sources in Chinese Spring for each of Lancer's genes")
my_genes <- cs_id_filler(data = my_genes, library = cs_id, chr = chromosome, ref.var = tolower(reference_cultivar))
print("2.7. Extract the names of the genes in chromosome 5B separated by commas") 
my_gene_sources <- as.character(my_genes$cs_id[!grepl("^source", my_genes$cs_id)])
my_gene_sources_text <- paste(my_gene_sources, collapse = ", ")
write.table(x = my_gene_sources_text, file = "my_gene_sources.txt", sep = "", row.names = F, col.names = F, quote = F)
my_variety_genes <- as.character(my_genes$var_id[!grepl("^source", my_genes$cs_id)])
my_variety_genes_text <- paste(my_variety_genes, collapse = ", ")
write.table(x = my_variety_genes_text, file = "my_variety_genes.txt", sep = "", row.names = F, col.names = F, quote = F)
print_my_genes <- data.frame( "var_id" = my_genes$var_id, "chinese_id" = my_genes$cs_id)
write.csv2(x = print_my_genes, file = "my_genes_and_sources.csv", row.names = F, quote = F)


############### 3. PROVE IF THE HAPLOBLOCK WAS CALLED BY GENE-BASED BLAST PAIRWISE ALIGNMENTS ###############   


print("3. PROVE IF THE HAPLOBLOCK WAS CALLED BY GENE-BASED BLAST PAIRWISE ALIGNMENTS")   
print("3.1. Download the zip file, decompress it in the working directory and read the tables") 
HC_gtf <- read.table("IWGSC_v1.2_HC_20200615.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
LC_gtf <- read.table("IWGSC_v1.2_LC_20200615.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
ALL_gtf <- rbind(HC_gtf, LC_gtf)
print("3.2. Extract the BLAST alignments and put them in table with Chinese Spring genes and their location in the IWGSC genome (long process)") 
BLAST_library <- read_pairwise_position(blast_path_gz = "varieties_all_identities_2000bp.tab.gz", gtf = ALL_gtf, write_table = "BLAST_library.tab")
print("3.3. Extract a subset with Lancer-Paragon comparison in chr5B") 
BLAST_subset <- BLAST_library[grepl(paste0(tolower(reference_cultivar), "->", tolower(query_cultivar), sep = ""), BLAST_library$aln_type) & grepl(chromosome, BLAST_library$chr),]
print("3.4. Use the vector with the genes identified in the previous step to extract another subset with only the genes in our haploblock") 
BLAST_subset <- BLAST_subset[grepl(paste(my_gene_sources, collapse = "|"), BLAST_subset$transcript),]
write.csv2(x = BLAST_subset, file = "BLAST_subset.csv", row.names = F, quote = F)
colnames(BLAST_subset)[1] <- "cs_transcript"
colnames(BLAST_subset)[14] <- "cs_start"
colnames(BLAST_subset)[15] <- "cs_end"
BLAST_subset$var_transcript <- my_genes$var_id[grepl(paste(BLAST_subset$cs_transcript, collapse = "|"), my_genes$cs_id)]
BLAST_subset$var_start <- my_genes$start[grepl(paste(BLAST_subset$var_transcript, collapse = "|"), my_genes$var_id)]
BLAST_subset$var_end <- my_genes$end[grepl(paste(BLAST_subset$var_transcript, collapse = "|"), my_genes$var_id)]
my_gene_sources[my_gene_sources %in% BLAST_subset$cs_transcript == FALSE]
my_genes$var_id[grepl(paste(my_gene_sources[my_gene_sources %in% BLAST_subset$cs_transcript == FALSE], collapse = "|"), my_genes$cs_id)]
my_gene_sources_filtered_by_Brinton <- my_gene_sources[my_gene_sources %in% BLAST_subset$cs_transcript == TRUE]
my_variety_genes_filtered_by_Brinton <- my_genes$var_id[grepl(paste(my_gene_sources [my_gene_sources %in% BLAST_subset$cs_transcript == TRUE], collapse = "|"), my_genes$cs_id)]
my_gene_sources_filtered_by_Brinton_text <- paste(my_gene_sources_filtered_by_Brinton, collapse = ", ")
write.table(x = my_gene_sources_filtered_by_Brinton_text, file = "my_gene_sources_filtered_by_Brinton_text.txt", sep = "", row.names = F, col.names = F, quote = F)
my_variety_genes_filtered_by_Brinton_text <- paste(my_variety_genes_filtered_by_Brinton, collapse = ", ")
write.table(x = my_variety_genes_filtered_by_Brinton_text, file = "my_variety_genes_filtered_by_Brinton_text.txt", sep = "", row.names = F, col.names = F, quote = F)
print_blast <- data.frame( "var_id" = my_variety_genes_filtered_by_Brinton_text, "chinese_id" = my_gene_sources_filtered_by_Brinton_text)
write.csv2(x = print_blast, file = "my_genes_and_sources_filtered_by_Brinton.csv", row.names = F, quote = F)
print("3.6. Calculate the percentage of identity in windows of 20 genes where genes containing Ns are filtered out") 
BLAST_subset <- BLAST_subset[ BLAST_subset$Ns_total == 0, ]
window_BLAST_subset <- edited_calculate_pid_windows(aln_data = BLAST_subset)
blocks_BLAST_subset <- assign_blocks(window_BLAST_subset)
print(blocks_BLAST_subset)