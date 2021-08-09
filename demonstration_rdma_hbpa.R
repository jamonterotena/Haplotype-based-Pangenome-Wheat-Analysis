########################################        HAPLOTYPE-BASED PANGENOME ANALYSIS IN WHEAT        ########################################         


# As the writer of this script I would like to express appreciation and give all credit to the authors of the following article for creating original raw data and scripts, which are required for this R script:

# Brinton, J., Ramirez-Gonzalez, R. H., Simmonds, J., Wingen, L., Orford, S., Griffiths, S., 10 Wheat Genome Project, Haberer, G., Spannagl, M., Walkowiak, S., Pozniak, C., & Uauy, C. (2020). A haplotype-led approach to increase the precision of wheat breeding. Communications biology, 3(1), 712. https://doi.org/10.1038/s42003-020-01413-2

# This script was written in R 4.1.0 by Jose Antonio Montero Tena as part of his master thesis in 2021. Free use and modification of this script for personal interest is encouraged.



# INTRODUCTION: SNP-haplotypes discovered from marker arrays can be redefined as sequence-based haplotypes, which are identical-by-state sequences sharing 100 percentage of identity between the genomes and discovered by sequence comparison. Brinton et al developed a method in 2020 that defined haploblocks as physical regions with >= 99.99 pident across fixed-size bins of 5-, 2.5- or 1-Mbp between pairwise comparisons of assemblies in the 10+ Wheat Genome Project. This method will subsequently be referred to as the Brinton's method. In these sequence-based haploblocks, haplotypes can be either shared haplotypes, when pairs of assemblies share one identical-by-state sequence in this region, or unique haplotypes, when no other assembly presents the same sequence in the haploblock. Additionally, Brinton and colleagues published crop-haplotypes.com, a website that provides an interactive graphic visualization of the shared haplotypes between the wheat genome assemblies. These assemblies were chromosome-level genome assemblies of nine wheat assemblys (ArinaLrFor, Jagger, Julius, Lancer, Landmark, Mace, Norin61, Stanley, SY-Mattis) and the Chinese Spring RefSev1.0 assembly alongside scaffold-level assemblies of five additional assemblys (Cadenza, Claire, Paragon, Robigus and Weebill).
# AIM: The aim of this script is to use Brinton's method and resources to conduct both a chromosome-scale and small-scale analysis of the pairwise comparisons between assemblies. The small-scale analysis is an original method to this script and allows the analysis of target regions within chromosomes. The aim of the small-scale analysis is to map physical start and end coordinates of the predicted haploblock both in the reference and query genomes by considering the region coverage with alignments, in other words, what is the number of alignments the haplotype prediction is based on and if there is any discontinuity that could cause false positive predictions. Also, this script aims to identify the genes contained in mapped haploblock, both in the haplotype-carrying genomes and in the IWGSC Chinese Spring Ref v1.1. Brinton et al called haplotypes both from mummer and gene-based BLAST pairwise alignments and selected for matching positions the longer predictions. Both types of analysis are provided in this script, although the results from these methods do not always match.
# BACKGROUND: The SNP-haplotype Hap-5B-RDMa-h2, associated with increased root dry mass (RDM), is taken as an example. This haplotype was first found by GWAS upon 9 SNP marker alleles located in the wheat chromosome 5B and matched with the assemblies of the assemblies LongReach Lancer and Paragon by BLAST. SNP-RDMa-h2 BLAST position in Lancer's chromosome 5B (655760000-656600000 bp) was observed in crop-haplotypes.com in order to tell if the SNP-haploblock region matched with the region of any sequence-based haploblock, discovered with Brinton's method, in this area. The following blocks were observed in the pairwise comparison Lancer-Paragon under each bin size (bin size - start - end - length): 5-Mbp - 660.000.000 - 665.000.000 - 5.000.000 // 2.5-Mbp - 657.500.000 - 662.500.000 - 5.000.000 // 1-Mbp - 656.000.000 - 663.000.000 - 7.000.000 (start and end according to Lancer's coordinate system as Lancer acts as the reference in this comparison). Considering that blocks and their positions can vary between bin sizes, this observation was sufficient evidence to confirm the match between the SNP- and a potential sequence-based haplotype. Despite of this, Brinton's method must be applied on this redefined region to understand if the prediction made on crop-haplotypes.com is reliable. Additionally, other research questions as the exact start and end position of the sequence-based haploblock or its genes can be answered in this script.


##### RUNNING THE SCRIPT ##### 

# The script is designed so that you only have to edit only the upcoming parameters. You can run the script by lines (recommended for first time) or simply pressing 'Source' in the upper right corner of this window. The running time is expected to be between 5 and 7 minutes. Please, be patient and only look for solutions if error messages interrupt the pipeline. If errors do not stop coming up, you can contact the script editor at the e-mail adress jmonterotena@gmail.com. 


##### DEFINING PARAMETERS #####

chromosome <- "5B" # Copy the same format
reference_assembly <- "Lancer" # Alignment coordinates will apply for this assembly's chromosome. Can NEVER be a scaffold-level assembly (Cadenza, Claire, Paragon, Robigus, Weebil). Follow the next naming guidelines: "arinalrfor", "jagger", "julius", "lancer", "landmark", "mace", "norin61", "stanley", "sy_mattis", "chinese" (works with capital letters)
query_assembly <- "Paragon" # This genome is aligned against the reference chromosome. Can be a scaffold-level assembly. Follow the next naming guidelines: "cadenza", "claire", "paragon", "robigus", "weebil" (works with capital letters)
zoom_start <- 655000000 # Start of the region highlighted in the small-scale analysis, in which haploblocks are sought for. ADVICE: use multiples of 5 so that the bins match with those on crop-haplotypes
zoom_end <- 665000000 # End of the region highlighted in the small-scale analysis, in which haploblocks are sought for. ADVICE: use multiples of 5 so that the bins match with those on crop-haplotypes
target_start <- 655700000 # If you have evidence of a haplotype within your zoom region that you want to compare with the sequence-based haplotypes obtained in this script, write its start coordinate in the reference chromosome. If not, simply write NA. The target region must be within the zoom region
target_end <- 656600000 # If you have evidence of a haplotype within your zoom region that you want to compare with the sequence-based haplotypes obtained in this script, write its end coordinate in the reference chromosome. If not, simply write NA. The target region must be within the zoom region

selected_start <- 655760000 # You do not need to wrorry about this yet. Complete section one and then come back to define the new coordinates of your haploblocks. Genes will be extracted from this position. If you have no interest in redefining your region, simply write 'target_start' or 'zoom_start', to keep with the previous coordinates
selected_end <- 662740000 # You do not need to wrorry about this yet. Complete section one and then come back to define the new coordinates of your haploblocks. Genes will be extracted until this position. If you have no interest in redefining your region, simply write 'target_end' or 'zoom_end', to keep with the previous coordinates

##### Running functions ##### 

# Source the script 'functions.hbpa.r' to read the packages and functions required for this script.

############### 1. CALLING HAPLOTYPES FROM RAW DATA CONTAINING >= 20-KBP-LONG MUMMER PAIRWISE ALIGNMENTS BETWEEN LANCER AND PARAGON ###############   

# Requirements:

# - Raw delta files in 'all_20_kb_filtered_delta_CHROMOSOME_tables.rds'; in this case the chromosome is 5B (downloadable from https://opendata.earlham.ac.uk/wheat/under_license/toronto/Brinton_etal_2020-05-20-Haplotypes-for-wheat-breeding/nucmer/). Brinton et al filtered out alignments under 20000Kbp of length in order to get rid of non-syntenic retrotransposons
# - Script 'functions.hbpa.r' (downloadable from https://github.com/MonteroJLU/Haplotype-based-Pangenome-Wheat-Analysis.git) 
# - Text file 'table_assembly_chr_length.txt' (downloadable from https://github.com/MonteroJLU/Haplotype-based-Pangenome-Wheat-Analysis.git)

##########  1.1. CHROMOSOME-SCALE ANALYSIS ##########  

##### 1.1.1. Download and save the required documents in the working directory ##### 

##### 1.1.2. Run the script 'functions_hbpa.r' ##### 

##### 1.1.3. Read the rds file into a data frame: ##### 

aln_library <- readRDS(file = paste0("all_20_kb_filtered_delta_", chromosome,"_tables.rds"))

##### 1.1.4. Create a subset of the data frame containing only the alignments from the pairwise comparison Lancer-Paragon, where Lancer acts as the the reference genome: ##### 

aln_subset <- aln_library[grepl(paste0("^", tolower(reference_assembly), sep = ""), aln_library$comparison) & grepl(tolower(query_assembly), aln_library$comparison),]

# Description of the column headers: rs: reference start, re: reference end, qs: query start, qe: query end, error: number of unmatches, qid: query identification, rid: reference identification, strand: forward or reverse strand, r_length: length of the alignment in the reference, perc_id: percentage of identity, perc_id_factor: factor that summarises perc_id, r_mid: midpoint the alignment in reference ((r_end-r_start)/2), q_mid: query midpoint, comparison: assemblys compared, chrom: chromosome.

##### 1.1.5. Scatter-plot the alignment midpoints (X: r_mid, Y: q_mid) ##### 

plot_diagonal_scatterplot(aln_subset, cap_lower = 90.00, cap_upper = 100, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 100000000)

# Notice unexpected vertical lines in the scatter-plot due to the use of a scaffold-level assembly as query (CLA-SLA comparison). However, this function works well for CLA-CLA comparisons.

##### 1.1.6. Dot-plot the alignments to show percentage of identity and alignment length (X: r_mid, Y: perc_id) ##### 

plot_aln_pid_and_length(aln_subset, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 100000000, dot_size = 2)

##### 1.1.7. Check for alignment properties in the data frame (percentage of alignment coverage of reference chromosome, average aln length and av aln number) ##### 

table_chrlength_v_assembly <- read.table(file = "table_assembly_chr_length.txt", sep = "\t", header = TRUE)
chr_length <- table_chrlength_v_assembly$sequence_length[grepl(tolower(reference_assembly), table_chrlength_v_assembly$assembly_name) & grepl(paste0("chr", unique(aln_subset$chrom)), table_chrlength_v_assembly$sequence_name)]
print(paste0("chr", unique(aln_subset$chrom), " is ", chr_length, " bp-long in ", reference_assembly))
print(paste0(round(mean(aln_subset$r_length), 0), " is the average alignment length in ", reference_assembly, "-", query_assembly, paste0(" chr", unique(aln_subset$chrom)), " comparison"))
print(paste0(nrow(aln_subset), " is the number of alignments in ", reference_assembly, "-", query_assembly, paste0(" chr", unique(aln_subset$chrom)), " comparison")) 
print(paste0(round((sum(aln_subset$r_length)/chr_length*100), 0), "% is the alignment coverage in ", reference_assembly, "-", query_assembly, paste0(" chr", unique(aln_subset$chrom)), " comparison")) 
coverage <- COV(aln_library, chr_length = table_chrlength_v_assembly)
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

##### 1.1.8. Boxplot the alignments (X: bin, Y: perc_id_median) ##### 

plot_boxplots_bin_median(aln_subset, bin_size = 10000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 50000000, show_outliers = FALSE)
plot_boxplots_bin_median(aln_subset, bin_size = 5000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 50000000, show_outliers = FALSE)
plot_boxplots_bin_median(aln_subset, bin_size = 2500000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 50000000, show_outliers = FALSE)
plot_boxplots_bin_median(aln_subset, bin_size = 1000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 50000000, show_outliers = FALSE)

##### 1.1.9. Plot the median percentage of identity across the chromosome or median line (X: r_mid, Y: perc_id_median) ##### 

medians_aln_subset_10Mbp <- plot_line_bin_median(aln_subset, bin_size = 10000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 50000000)
medians_aln_subset_10Mbp
medians_aln_subset_5Mbp <- plot_line_bin_median(aln_subset, bin_size = 5000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 50000000)
medians_aln_subset_5Mbp
medians_aln_subset_2.5Mbp <- plot_line_bin_median(aln_subset, bin_size = 2500000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 50000000)
medians_aln_subset_2.5Mbp
medians_aln_subset_1Mbp <- plot_line_bin_median(aln_subset, bin_size = 1000000, bin_start = 0, bin_end = max(aln_subset$re), cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 50000000)
medians_aln_subset_1Mbp

##### 1.1.10. Print information about the haploblock predictions in the pairwise comparison across the reference chromosome ##### 

medians_aln_subset_10Mbp_bin_info <- assign_blocks_mummer(median_cutoffs = medians_aln_subset_10Mbp[["data"]], original_file = aln_subset)
medians_aln_subset_5Mbp_bin_info <- assign_blocks_mummer(median_cutoffs = medians_aln_subset_5Mbp[["data"]], original_file = aln_subset)
medians_aln_subset_2.5Mbp_bin_info <- assign_blocks_mummer(median_cutoffs = medians_aln_subset_2.5Mbp[["data"]], original_file = aln_subset)
medians_aln_subset_1Mbp_bin_info <- assign_blocks_mummer(median_cutoffs = medians_aln_subset_1Mbp[["data"]], original_file = aln_subset)
medians_aln_subset_10Mbp_block_summary <- block_summary(medians_aln_subset_10Mbp_bin_info, bin_size = 10000000, reference_name = reference_assembly, query_name = query_assembly )
medians_aln_subset_5Mbp_block_summary <- block_summary(medians_aln_subset_5Mbp_bin_info, bin_size = 5000000, reference_name = reference_assembly, query_name = query_assembly )
medians_aln_subset_2.5Mbp_block_summary <- block_summary(medians_aln_subset_2.5Mbp_bin_info, bin_size = 2500000, reference_name = reference_assembly, query_name = query_assembly )
medians_aln_subset_1Mbp_block_summary <- block_summary(medians_aln_subset_1Mbp_bin_info, bin_size = 1000000, reference_name = reference_assembly, query_name = query_assembly )
print(medians_aln_subset_10Mbp_block_summary)
print(medians_aln_subset_5Mbp_block_summary)
print(medians_aln_subset_2.5Mbp_block_summary)
print(medians_aln_subset_1Mbp_block_summary)

##########  1.2. SMALL-SCALE ANALYSIS ##########  

##### 1.2.1. Scatter-plot the alignment midpoints across the zoom region (X: r_mid, Y: q_mid) ##### 

plot_diagonal_scatterplot(aln_subset, xmin = zoom_start, xmax = zoom_end, cap_lower = 90.00, cap_upper = 100, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 1000000)

##### 1.2.2. Dot-plot the alignments to show percentage of identity and alignment length in the zoom region (X: r_mid, Y: perc_id) ##### 

graph <- plot_aln_pid_and_length(data = aln_subset, xmin = zoom_start, xmax = zoom_end, ymin = 98.5, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 1000000, dot_size = 4)
graph

##### 1.2.3. Check for alignment properties in the zoom region #####

aln_target <- aln_subset[(aln_subset$r_mid >= target_start) & (aln_subset$r_mid <= target_end),]
print(paste0(round(mean(aln_target$r_length), 0), " is the average alignment length for the target region between ", target_start/1e06, " and ", target_end/1e06, " Mbp in ", reference_assembly, "-", query_assembly, paste0(" chr", unique(aln_target$chrom)), " comparison")) 
print(paste0(nrow(aln_target), " is the number of alignments for the target region between ", target_start/1e06, " and ", target_end/1e06, " Mbp in ", reference_assembly, "-", query_assembly, paste0(" chr", unique(aln_target$chrom)), " comparison"))
print(paste0(round((sum(aln_target$r_length)/(target_end-target_start)*100), 0), "% is the alignment coverage for the target region between ", target_start/1e06, " and ", target_end/1e06, " Mbp in ", reference_assembly, "-", query_assembly, paste0(" chr", unique(aln_target$chrom)), " comparison") )
aln_zoom <- aln_subset[(aln_subset$r_mid >= zoom_start) & (aln_subset$r_mid <= zoom_end),]
print(paste0(round(mean(aln_zoom$r_length), 0), " is the average alignment length for the zoom region between ", zoom_start/1e06, " and ", zoom_end/1e06, " Mbp in ", reference_assembly, "-", query_assembly, paste0(" chr", unique(aln_subset$chrom)), " comparison")) 
print(paste0(nrow(aln_zoom), " is the number of alignments for the zoom region between ", zoom_start/1e06, " and ", zoom_end/1e06, " Mbp in ", reference_assembly, "-", query_assembly, paste0(" chr", unique(aln_subset$chrom)), " comparison"))
print(paste0(round((sum(aln_zoom$r_length)/(zoom_end-zoom_start)*100), 0), "% is the alignment coverage for the zoom region between ", zoom_start/1e06, " and ", zoom_end/1e06, " Mbp in ", reference_assembly, "-", query_assembly, paste0(" chr", unique(aln_subset$chrom)), " comparison")) 
bin_size <- c(5000000, 2500000, 1000000)
names(bin_size) <- c("bin size: 5-Mbp", "bin size: 2.5-Mbp", "bin size: 1-Mbp")
for (i in 1:3){
  print("average expected number of alignments per bin across zoom region")
  print(nrow(aln_zoom)/((zoom_end-zoom_start)/bin_size[i]))
}

##### 1.2.4. Boxplot the bin median percentage of identity in the zoom region to check for outliers #####

plot_boxplots_bin_median(aln_subset, bin_size = 5000000, bin_start = zoom_start, bin_end = zoom_end, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 1000000, show_outliers = TRUE)
plot_boxplots_bin_median(aln_subset, bin_size = 2500000, bin_start = zoom_start, bin_end = zoom_end, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 1000000, show_outliers = TRUE)
plot_boxplots_bin_median(aln_subset, bin_size = 1000000, bin_start = zoom_start, bin_end = zoom_end, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 1000000, show_outliers = TRUE)

##### 1.2.5. Plot the median line in the zoom region to see the haploblock predictions at different bin sizes #####

plot_line_bin_median(aln_subset, bin_size = 5000000, bin_start = zoom_start, bin_end = zoom_end, ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 1000000)
plot_line_bin_median(aln_subset, bin_size = 2500000, bin_start = zoom_start, bin_end = zoom_end, ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 1000000)
plot_line_bin_median(aln_subset, bin_size = 1000000, bin_start = zoom_start, bin_end = zoom_end, ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 1000000)

##### 1.2.6. Dot-plot the zoom region with the haploblock predictions and the target region and make decisions regarding the start and end limits of the sequence-based haploblock ##### 

target <- data.frame(target_start, target_end)
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 5000000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 2500000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 1000000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)

##### DEFINING NEW PARAMETERS ##### 

selected_start <- 655760000 # Genes will be extracted from this position. If you have no interest in redefining your region, simply write 'target_start' or 'zoom_start', to keep with the previous coordinates
selected_end <- 662740000 # Genes will be extracted until this position. If you have no interest in redefining your region, simply write 'target_end' or 'zoom_end', to keep with the previous coordinates
target_text <- "Potential new haploblock" # Text to print on the selected region

selected_haploblock <- data.frame(selected_start, selected_end)
print(selected_haploblock)
##### 1.2.7. Plot summary graphs ##### 

# This would be the physical location and length of the sequence-based haploblock Hap-5B-RDMa-h2 - start: 655.760.000 bp - end: 662.740.000 bp - length: 6.980.000 bp (reference: Lancer, query: Paragon). This information is shown in the following summary graphs, including the names of Paragon's scaffolds cut by the limits:

plot_aln_and_bins(print_tables = FALSE, aln_subset = aln_subset, bin_size = 1000000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = selected_haploblock, target_text = target_text, fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
plot_bins_and_selected_region(print_tables = FALSE, aln_subset = aln_subset, bin_size = 1000000, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = selected_haploblock, target_text = target_text, fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
plot_aln_and_bins(print_tables = FALSE, aln_subset = aln_subset, bin_size = 5000000, zoom_start = 0, zoom_end = max(aln_subset$re), highlighted_target = selected_haploblock, target_text = target_text, fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 100000000, prediction_text = FALSE, aln_text = F)
plot_bins_and_selected_region(print_tables = FALSE, aln_subset = aln_subset, bin_size = 5000000, zoom_start = 0, zoom_end = max(aln_subset$re), highlighted_target = selected_haploblock, target_text = target_text, fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 100000000, prediction_text = FALSE)



############### 2. IDENTIFY GENES WITHIN THE JUST-MAPPED HAPLOBLOCK ###############   

# Requirements:

# - File 'projectedGenes__Triticum_aestivum_REFERENCEassembly_v1.0.gff', in this case the reference is LongReach Lancer (downloadable from https://webblast.ipk-gatersleben.de/downloads/wheat/gene_projection/)
# - File 'geneid_2_chinese.sourceid.txt' (downloadable from https://webblast.ipk-gatersleben.de/downloads/wheat/gene_projection/)

##### 2.1. Download the files and save them in the working directory #####

##### 2.2. Read the gff file containing the gene model projection for the reference assembly #####

ref_gff <- read.table ( file = "projectedGenes__Triticum_aestivum_LongReach_Lancer_v1.0.gff", sep = "\t" , header = F, stringsAsFactors = F)

##### 2.3. Create subsets for gene projections only and edit the table #####

ref_gff_gene_only <- ref_gff[ref_gff$V3 == "gene",]
colnames(ref_gff_gene_only) <- c("chr", "annotation", "biotype", "start", "end", "score", "strand", "info", "ref_id")
ref_gff_gene_only$var_id <- gsub("ID=", "", ref_gff_gene_only$ref_id)

##### 2.4. Create a subset for the haploblock #####

my_genes <- ref_gff_gene_only[grep(chromosome, ref_gff_gene_only$chr),]
my_genes <- my_genes[(my_genes$start >= selected_start) & (my_genes$end <= selected_end),]

##### 2.5. Download 'geneid_2_chinese.sourceid.txt' from the last link and read the table ##### 

cs_id <- read.table(file = "geneid_2_chinese.sourceid.txt", sep="\t", header=T, stringsAsFactors = F)

# This text file matches the gene projection with their sources in Chinese Spring annotation

##### 2.6. Add a column to the subset with the IDs of the gene sources in Chinese Spring for each of Lancer's genes ##### 

my_genes <- cs_id_filler(data = my_genes, library = cs_id, chr = chromosome, ref.var = tolower(reference_assembly))

##### 2.7. Extract the names of the genes in chromosome 5B separated by commas ##### 

my_gene_sources <- as.character(my_genes$cs_id[!grepl("^source", my_genes$cs_id)])
my_gene_sources_text <- paste(my_gene_sources, collapse = ", ")
write.table(x = my_gene_sources_text, file = "my_gene_sources.txt", sep = "", row.names = F, col.names = F, quote = F)

my_variety_genes <- as.character(my_genes$var_id[!grepl("^source", my_genes$cs_id)])
my_variety_genes_text <- paste(my_variety_genes, collapse = ", ")
write.table(x = my_variety_genes_text, file = "my_variety_genes.txt", sep = "", row.names = F, col.names = F, quote = F)

print_my_genes <- data.frame( "var_id" = my_genes$var_id, "chinese_id" = my_genes$cs_id)
write.csv2(x = print_my_genes, file = "my_genes_and_sources.csv", row.names = F, quote = F)



############### 3. PROVE IF THE HAPLOBLOCK WAS CALLED BY GENE-BASED BLAST PAIRWISE ALIGNMENTS ###############   

# Requirements:

# - File 'varieties_all_identities_2000bp.tar.gz' (downloadable from https://opendata.earlham.ac.uk/wheat/under_license/toronto/Brinton_etal_2020-05-20-Haplotypes-for-wheat-breeding/pairwise_blast/)
# - File 'iwgsc_refseq_v1.2_gene_annotation.zip' (downloadable from https://urgi.versailles.inrae.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.2/)

##### 3.1. Download the zip file, decompress it in the working directory and read the tables ##### 

HC_gtf <- read.table("IWGSC_v1.2_HC_20200615.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
LC_gtf <- read.table("IWGSC_v1.2_LC_20200615.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
ALL_gtf <- rbind(HC_gtf, LC_gtf)

##### 3.2. Extract the BLAST alignments and put them in table with Chinese Spring genes and their location in the IWGSC genome (long process) ##### 

BLAST_library <- read_pairwise_position(blast_path_gz = "varieties_all_identities_2000bp.tab.gz", gtf = ALL_gtf, write_table = "BLAST_library.tab")

##### 3.3. Extract a subset with Lancer-Paragon comparison in chr5B ##### 

BLAST_subset <- BLAST_library[grepl(paste0(tolower(reference_assembly), "->", tolower(query_assembly), sep = ""), BLAST_library$aln_type) & grepl(chromosome, BLAST_library$chr),]

##### 3.4. Use the vector with the genes identified in the previous step to extract another subset with only the genes in our haploblock ##### 

BLAST_subset <- BLAST_subset[grepl(paste(my_gene_sources, collapse = "|"), BLAST_subset$transcript),]
write.csv2(x = BLAST_subset, file = "BLAST_subset.csv", row.names = F, quote = F)

# Brinton et al (2020) only retained only gene projections consistent with the expected chromosome. This explains why this subset contains less genes that the total amount of annotated genes in our region

# Add the coordinates and transcripts from the variety #

colnames(BLAST_subset)[1] <- "cs_transcript"
colnames(BLAST_subset)[14] <- "cs_start"
colnames(BLAST_subset)[15] <- "cs_end"

BLAST_subset$var_transcript <- my_genes$var_id[grepl(paste(BLAST_subset$cs_transcript, collapse = "|"), my_genes$cs_id)]
BLAST_subset$var_start <- my_genes$start[grepl(paste(BLAST_subset$var_transcript, collapse = "|"), my_genes$var_id)]
BLAST_subset$var_end <- my_genes$end[grepl(paste(BLAST_subset$var_transcript, collapse = "|"), my_genes$var_id)]

# Notice that some genes were filtered out if they contained more than one projection in the expected chromosome, so the amount of genes shown in this table is expected to be smaller than those found in lancer's projected genes in the region. The genes that were excluded can be seen here:

my_gene_sources[my_gene_sources %in% BLAST_subset$cs_transcript == FALSE]
my_genes$var_id[grepl(paste(my_gene_sources[my_gene_sources %in% BLAST_subset$cs_transcript == FALSE], collapse = "|"), my_genes$cs_id)]

# We can extract another list containing only the genes that were present among the BLAST alignments

my_gene_sources_filtered_by_Brinton <- my_gene_sources[my_gene_sources %in% BLAST_subset$cs_transcript == TRUE]
my_variety_genes_filtered_by_Brinton <- my_genes$var_id[grepl(paste(my_gene_sources [my_gene_sources %in% BLAST_subset$cs_transcript == TRUE], collapse = "|"), my_genes$cs_id)]
my_gene_sources_filtered_by_Brinton_text <- paste(my_gene_sources_filtered_by_Brinton, collapse = ", ")
write.table(x = my_gene_sources_filtered_by_Brinton_text, file = "my_gene_sources_filtered_by_Brinton_text.txt", sep = "", row.names = F, col.names = F, quote = F)

my_variety_genes_filtered_by_Brinton_text <- paste(my_variety_genes_filtered_by_Brinton, collapse = ", ")
write.table(x = my_variety_genes_filtered_by_Brinton_text, file = "my_variety_genes_filtered_by_Brinton_text.txt", sep = "", row.names = F, col.names = F, quote = F)

print_blast <- data.frame( "var_id" = my_variety_genes_filtered_by_Brinton_text, "chinese_id" = my_gene_sources_filtered_by_Brinton_text)
write.csv2(x = print_blast, file = "my_genes_and_sources_filtered_by_Brinton.csv", row.names = F, quote = F)

##### 3.5. Calculate the percentage of identity in windows of 20 genes where genes containing Ns are filtered out ##### 

BLAST_subset <- BLAST_subset[ BLAST_subset$Ns_total == 0, ]

window_BLAST_subset <- edited_calculate_pid_windows(aln_data = BLAST_subset)
blocks_BLAST_subset <- assign_blocks(window_BLAST_subset)
print(blocks_BLAST_subset)

# We can see how in a region where a haploblock was assigned from mummer pairwise alignments at different bin sizes (5, 2.5, 1Mbp), no block was assigned when using gene-based BLAST pairwise alignments +/- 2000 bp flanking sequence. Originally, the blocks were called when they were confirmed by any of the methods and the longest prediction was taken, so this explains why crop-haplotypes.com shows haplotypes in this region. These two methods to assign blocks from MUMmer or BLAST pairwise alignments utilize different input and criteria that could result in differences regarding the blocks called:
# 1. The first receives as input any alignment between the reference and the query, however, the second method only utilizes alignments from gene model projections in the reference and query with 2000 bp up and downstream region
# 2. The first method filters out alignments under 20000 bp, meanwhile the second works with BLAST alignments with length 4000 pb + gene, so  the first method could technically filter out the alignments that the second method would use if they were not part of longer alignments
# 3. The second method filters out alignments with N's, which is not followed by the first method. Instead, the first method uses a lower
# cutoff 99.99% to accomodate for the presence of N's in the alignments
# 4. The cutoff for the first method is 99.99%, while the second method 100%, although drops out the two alignments with the lowest pidents by window. If the cutoff for the second method was the same as for the first method, blocks would have been called

# As a conclusion, the second method is very sensitive about sequencing quality and this seems to explaining how no blocks were called.