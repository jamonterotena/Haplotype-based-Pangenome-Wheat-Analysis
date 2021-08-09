###############   REQUIRED PACKAGES ###############   

library(plyr)
library(dplyr)
library(magrittr)
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(viridis)
library(stringr)
library("ggdendro")
library(reshape2)
library("grid")

###############   FUNCTIONS ###############   

##########  1. CALLING HAPLOTYPES FROM RAW DATA CONTAINING >= 20-KBP-LONG MUMMER PAIRWISE ALIGNMENTS BETWEEN LANCER AND PARAGON ##########  

readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

edited_bin_data <- function(data, bin_size, bin_start = 0, bin_end){
  bins <- seq(bin_start, bin_end, by = bin_size)
  data$bin <- NA
  for (i in bins){
    data$bin <- ifelse(((data$r_mid > (i-bin_size)) & (data$r_mid < (i-1))), i, data$bin)
  }
  return(data)
}

plot_diagonal_scatterplot <- function(data, xmin = 0, xmax = max(data$re), cap_lower, cap_upper, reference_name = data$rid, query_name = data$qid, x_label_gap){
  ggplot(data[data$rs > xmin & data$re < xmax,], aes(x=rs, xend=re, y=qs, yend=qe, colour=perc_id, shape = strand)) + geom_segment() +
    ggtitle(paste0("Reference vs query (diagonal scatterplot): ", reference_name, " vs ", query_name, ", zoom region: ", xmin/1000000, "-", xmax/1000000, " Mbp")) +
    geom_point(alpha=.5) +
    theme_bw() + 
    theme(strip.text.y=element_text(angle=180, size=5), strip.background=element_blank()) +
    xlab(paste0 ("alignment position in ", reference_name, "'s chr", unique(data$chrom), " (bp)")) + ylab(paste0("alignment position in ", query_name, "'s chr", unique(data$chrom), " (bp)")) +
    scale_x_continuous(breaks = seq( xmin, xmax, by = x_label_gap), position = "bottom") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_viridis(limits = c(cap_lower, cap_upper))
}

plot_aln_pid_and_length <- function(data, xmin = 0, xmax = max(data$re), ymin = 97, ymax = 100, reference_name = data$rid, query_name = data$qid, x_label_gap = 5000000, dot_size = 2){
  ggplot(data[data$r_mid > xmin & data$r_mid < xmax,], aes(x=r_mid, y=perc_id, colour=r_length)) +
    theme_bw() + 
    xlab(paste0 ("alignment position in ", reference_name, "'s chr", unique(data$chrom), " (bp)")) + 
    ylab(paste0('alignment percentage of identity vs ', query_name)) +
    geom_point(alpha = .5, size = dot_size) +
    ylim(ymin, ymax+0.05) +
    scale_colour_viridis() +
    geom_hline(yintercept = 99.99, linetype = "dashed", color = "red") +
    annotate ("text", colour = "red", size = 4.5, x = xmin, y = 100.025, label = "cutoff 99.99%" ) +
    scale_x_continuous(breaks = seq( xmin, xmax, by = x_label_gap), position = "top") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0("Percentage of identity of individual alignments: ", " ", reference_name, " vs ", query_name, ", zoom region: ", xmin/1e06, "-", xmax/1e06, " Mbp"))
}

COV <- function(data = data.frame(), chr_length = table_assembly_chr_length, assembly_names = c("arinalrfor", "jagger", "julius", "lancer", "landmark", "mace", "norin61", "stanley", "sy_mattis", "chinese"), scaffold_level_assemblies = c("claire", "cadenza", "paragon", "robigus", "weebil")){
  all_scaffold_level_assemblies = c("claire", "cadenza", "paragon", "robigus", "weebil")
  cov_assembly_CLA_CLA <- vector()
  cov_assembly_CLA_SLA <- vector()
  for (i in 1:length(assembly_names)){
    to_add_CLA_CLA <- sum(data$r_length[grepl(paste("^", assembly_names[i], sep = ""), data$comparison) & !grepl(paste(all_scaffold_level_assemblies, collapse = "|"), data$comparison) & !grepl("spelta", data$comparison)])/(length(assembly_names)-1)
    to_add_CLA_SLA <- sum(data$r_length[grepl(paste("^", assembly_names[i], sep = ""), data$comparison) & grepl(paste(scaffold_level_assemblies, collapse = "|"), data$comparison) & !grepl("spelta", data$comparison)])/length(scaffold_level_assemblies)
    to_add_CLA_CLA <- to_add_CLA_CLA/chr_length$sequence_length[grepl(paste0("chr", data$chrom[1], collapse = ""), chr_length$sequence_name) & grepl(assembly_names[i], chr_length$assembly_name)]
    to_add_CLA_SLA <- to_add_CLA_SLA/chr_length$sequence_length[grepl(paste0("chr", data$chrom[1], collapse = ""), chr_length$sequence_name) & grepl(assembly_names[i], chr_length$assembly_name)]
    cov_assembly_CLA_CLA <- append(cov_assembly_CLA_CLA, to_add_CLA_CLA*100)
    cov_assembly_CLA_SLA <- append(cov_assembly_CLA_SLA, to_add_CLA_SLA*100)
  }
  names(cov_assembly_CLA_CLA) <- assembly_names
  names(cov_assembly_CLA_SLA) <- assembly_names
  mean_cov_assembly_CLA_CLA <- mean(cov_assembly_CLA_CLA)
  mean_cov_assembly_CLA_SLA <- mean(cov_assembly_CLA_SLA)
  output <- list('Average % coverage with CHROMOSOME-LEVEL ASSEMBLIES as query' = mean_cov_assembly_CLA_CLA,
                 '(By assembly)' = cov_assembly_CLA_CLA,
                 'Average % coverage with SCAFFOLD-LEVEL ASSEMBLIES as query' = mean_cov_assembly_CLA_SLA,
                 '(By assembly)' = cov_assembly_CLA_SLA)
  return(output)
}

LENGTH <- function(data = data.frame(), assembly_names = c("arinalrfor", "jagger", "julius", "lancer", "landmark", "mace", "norin61", "stanley", "sy_mattis", "chinese"), scaffold_level_assemblies = c("claire", "cadenza", "paragon", "robigus", "weebil")){
  all_scaffold_level_assemblies = c("claire", "cadenza", "paragon", "robigus", "weebil")
  length_assembly_CLA_CLA <- vector()
  length_assembly_CLA_SLA <- vector()
  for (i in 1:length(assembly_names)){
    to_add_CLA_CLA <- mean(data$r_length[grepl(paste("^", assembly_names[i], sep = ""), data$comparison) & !grepl(paste(all_scaffold_level_assemblies, collapse = "|"), data$comparison) & !grepl("spelta", data$comparison)])
    to_add_CLA_SLA <- mean(data$r_length[grepl(paste("^", assembly_names[i], sep = ""), data$comparison) & grepl(paste(scaffold_level_assemblies, collapse = "|"), data$comparison) & !grepl("spelta", data$comparison)])
    length_assembly_CLA_CLA <- append(length_assembly_CLA_CLA, to_add_CLA_CLA)
    length_assembly_CLA_SLA <- append(length_assembly_CLA_SLA, to_add_CLA_SLA)
  }
  names(length_assembly_CLA_CLA) <- assembly_names
  names(length_assembly_CLA_SLA) <- assembly_names
  mean_length_assembly_CLA_CLA <- mean(length_assembly_CLA_CLA)
  mean_length_assembly_CLA_SLA <- mean(length_assembly_CLA_SLA)
  output <- list('Average aln length with CHROMOSOME-LEVEL ASSEMBLIES as query' = mean_length_assembly_CLA_CLA,
                 '(By assembly)' = length_assembly_CLA_CLA,
                 'Average aln length with SCAFFOLD-LEVEL ASSEMBLIES as query' = mean_length_assembly_CLA_SLA,
                 '(By assembly)' = length_assembly_CLA_SLA)
  return(output)
}

NUMBER <- function(data = data.frame(), assembly_names = c("arinalrfor", "jagger", "julius", "lancer", "landmark", "mace", "norin61", "stanley", "sy_mattis", "chinese"), scaffold_level_assemblies = c("claire", "cadenza", "paragon", "robigus", "weebil")){
  all_scaffold_level_assemblies = c("claire", "cadenza", "paragon", "robigus", "weebil")
  number_assembly_CLA_CLA <- vector()
  number_assembly_CLA_SLA <- vector()
  for (i in 1:length(assembly_names)){
    to_add_CLA_CLA <- nrow(data[grepl(paste("^", assembly_names[i], sep = ""), data$comparison) & !grepl(paste(all_scaffold_level_assemblies, collapse = "|"), data$comparison) & !grepl("spelta", data$comparison),])/(length(assembly_names)-1)
    to_add_CLA_SLA <- nrow(data[grepl(paste("^", assembly_names[i], sep = ""), data$comparison) & grepl(paste(scaffold_level_assemblies, collapse = "|"), data$comparison) & !grepl("spelta", data$comparison),])/length(scaffold_level_assemblies)
    number_assembly_CLA_CLA <- append(number_assembly_CLA_CLA, to_add_CLA_CLA)
    number_assembly_CLA_SLA <- append(number_assembly_CLA_SLA, to_add_CLA_SLA)
  }
  names(number_assembly_CLA_CLA) <- assembly_names
  names(number_assembly_CLA_SLA) <- assembly_names
  mean_number_assembly_CLA_CLA <- mean(number_assembly_CLA_CLA)
  mean_number_assembly_CLA_SLA <- mean(number_assembly_CLA_SLA)
  output <- list('Average aln number with CHROMOSOME-LEVEL ASSEMBLIES as query' = mean_number_assembly_CLA_CLA,
                 '(By assembly)' = number_assembly_CLA_CLA,
                 'Average aln number with SCAFFOLD-LEVEL ASSEMBLIES as query' = mean_number_assembly_CLA_SLA,
                 '(By assembly)' = number_assembly_CLA_SLA)
  return(output)
}

plot_boxplots_bin_median <- function(data, bin_size = 10000000, bin_start = 0, bin_end, cut_off = 99.99, ymin = 97, ymax = 100, reference_name = data$rid, query_name = data$qid, x_label_gap, show_outliers){
  comparison_filt_bin <- edited_bin_data(data, bin_size, bin_start, bin_end)
  comparison_medians <- aggregate(perc_id ~ bin, data = comparison_filt_bin, FUN=median)
  colnames(comparison_medians)[2] <- "perc_id_median"
  comparison_medians$cut_off <- NA
  comparison_medians$cut_off <- ifelse(comparison_medians$perc_id >= cut_off, paste0(">=",cut_off), paste0("<",cut_off))
  comparison_medians$cut_off <- as.factor(comparison_medians$cut_off)
  
  comparison_to_plot <- merge(comparison_filt_bin, comparison_medians)
  
  if (show_outliers == TRUE){
    ggplot(comparison_to_plot, aes(x=bin, y = perc_id, group = bin, fill = cut_off)) +
      geom_boxplot(outlier.shape = NA) +
      geom_dotplot(binaxis = "y", dotsize = 0.1, stackratio = 1, stackdir = "center", colour = "blue", binwidth = 0.1) + 
      ylim(ymin,ymax) +
      scale_fill_manual(values = c("73D055FF", "440154FF")) +
      labs(fill = "% id") +
      scale_colour_viridis() +
      xlab(paste0 ("bin position in ", reference_name, "'s chr", unique(data$chrom), " (bp)")) +
      ylab("Median bin percentage of identity") +
      scale_x_continuous(breaks = round( seq (bin_start, bin_end, by = x_label_gap), 1), position = "top") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle(paste0("Boxplot: ",reference_name, " vs ", query_name, ", zoom region: ", bin_start, "-", bin_end, " Mbp, bin size: ", bin_size/1000000, "-Mbp", ", cutoff: ", cut_off))
  } else {  ggplot(comparison_to_plot, aes(x=bin, y = perc_id, group = bin, fill = cut_off)) +
      geom_boxplot(outlier.shape = NA) +
      ylim(ymin,ymax) +
      scale_fill_manual(values = c("73D055FF", "440154FF")) +
      labs(fill = "% id") +
      scale_colour_viridis() +
      xlab(paste0 ("bin position in ", reference_name, "'s chr", unique(data$chrom), " (bp)")) +
      ylab("Median bin percentage of identity") +
      scale_x_continuous(breaks = round( seq (bin_start, bin_end, by = x_label_gap), 1), position = "top") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle(paste0("Boxplot: ",reference_name, " vs ", query_name, ", zoom region: ", bin_start, "-", bin_end, " Mbp, bin size: ", bin_size/1000000, "-Mbp", ", cutoff: ", cut_off))
  }
}

plot_line_bin_median <- function(data, bin_start = 0, bin_end, bin_size = 1e+07, cut_off = 99.99, ymin = 97, ymax = 100, reference_name = data$rid, query_name = data$qid, x_label_gap){
  comparison_filt_bin <- edited_bin_data(data = data, bin_size = bin_size, bin_start = bin_start, bin_end = bin_end)
  comparison_medians <- aggregate(perc_id ~ bin, data = comparison_filt_bin, FUN=median)
  colnames(comparison_medians)[2] <- "perc_id_median"
  comparison_medians$cut_off <- NA
  comparison_medians$cut_off <- ifelse(comparison_medians$perc_id >= cut_off, paste0(">=",cut_off), paste0("<",cut_off))
  comparison_medians$cut_off <- as.factor(comparison_medians$cut_off)
  
  ggplot(comparison_medians, aes(x=bin, y = perc_id_median, colour = cut_off)) +
    geom_line(colour = "grey") + 
    geom_point(size = 1)  + 
    ylim(ymin,ymax) +
    scale_colour_manual(values = c("73D055FF", "440154FF")) +
    labs(colour = "% id cutoff") +
    # scale_colour_viridis() +
    xlab(paste0 ("alignment position in ", reference_name, "'s chr", unique(data$chrom), " (bp)")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Median bin percentage of identity") +
    scale_x_continuous(breaks = round( seq (min(bin_start), max(bin_end), by = x_label_gap), 1), position = "top") +
    
    ggtitle(paste0("Median line: ", reference_name, " vs ", query_name, ", zoom region: ", bin_start, "-", bin_end, ", bin size: ", bin_size/1000000, "-Mbp"))
}

assign_blocks_mummer <-function(median_cutoffs, original_file){
  median_cutoffs_copy <- median_cutoffs
  median_cutoffs_copy$block_no <- NA
  block_no = 1
  
  for (i in seq(1, nrow(median_cutoffs_copy))){
    if(median_cutoffs_copy[i, "perc_id_median"] < 99.99){
      median_cutoffs_copy[i, "block_no"] <- "NO_BLOCK"
    } else if (median_cutoffs_copy[i, "perc_id_median"] >= 99.99){
      median_cutoffs_copy[i, "block_no"] <- block_no
      if (i > (nrow(median_cutoffs_copy)-3)){
      } else if ((median_cutoffs_copy[i+1, "perc_id_median"] < 99.99) & (median_cutoffs_copy[i+2, "perc_id_median"] < 99.99) & (median_cutoffs_copy[i+3, "perc_id_median"] < 99.99)){
        block_no <- block_no + 1
      }
    }
  }
  binSize <- (median_cutoffs_copy$bin[2] - median_cutoffs_copy$bin[1])
  median_cutoffs_copy$bin_start <- (median_cutoffs_copy$bin - binSize)
  median_cutoffs_copy$bin_end <- (median_cutoffs_copy$bin)
  number_block <- numeric()
  for (i in 1:nrow(median_cutoffs_copy)){
    tempo <- original_file[original_file$r_mid>=as.numeric(median_cutoffs_copy$bin_start[i]) & original_file$r_mid<=as.numeric(median_cutoffs_copy$bin_end[i]), ]
    ntempo <- nrow(tempo)
    number_block <- append(number_block, ntempo)
  }
  median_cutoffs_copy$aln_number <- number_block
  # print(paste0("BINS AT ", bin_size, "-MBP BIN SIZE"))
  # median_cutoffs$bin_size <- rep(binSize, nrow(median_cutoffs))
  return(median_cutoffs_copy)
}

block_summary <- function(median_cutoffs_copy, bin_size, reference_name = "NA", query_name ="NA", show_only_coords = FALSE){
  blocks <- unique(median_cutoffs_copy$block_no[!grepl("NO_BLOCK", median_cutoffs_copy$block_no)])
  blocks <- blocks[complete.cases(blocks)]
  
  block_positions <- data.frame(bin_size = character(), comparison = character(), block_no = numeric(), block_start = numeric(), block_end = numeric())
  
  for(block in blocks){
    bin_size_rep <- paste0(bin_size/1000000, "-Mbp", sep = "")
    block_data <- subset(median_cutoffs_copy, block_no == block)
    block_start <- (min(block_data$bin) - bin_size)
    block_end <- (max(block_data$bin))
    comparison <- paste0(reference_name, "->", query_name)
    to_add <- data.frame(bin_size = bin_size_rep,comparison = comparison, block_no = block, block_start = block_start, block_end = block_end)
    block_positions <- rbind(block_positions, to_add)
  }
  coords <- data.frame("start" = block_positions$block_start, "end" = block_positions$block_end)
  #print(paste0("BLOCK SUMMARY AT ", bin_size, "-MBP BIN SIZE"))
  ifelse(show_only_coords==TRUE, return(coords), return(block_positions))
}

plot_aln_and_bins <- function(aln_subset = data.frame(), bin_size = bin_size, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "target region", fill_target = "orange", color_target_text = "black", fill_predictions = "green", color_prediction_text = "black",  ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000, print_tables = TRUE, prediction_text = TRUE, aln_text = TRUE){
  medians_zoom <- plot_line_bin_median(data = aln_subset, bin_size = bin_size, bin_start = zoom_start, bin_end = zoom_end, ymin = ymin, cut_off = cut_off, reference_name = reference_assembly, query_name = query_assembly, x_label_gap = x_label_gap)
  medians_zoom_bin_info <- assign_blocks_mummer(medians_zoom[["data"]], original_file = aln_subset)
  aln_number_per_bin <- medians_zoom_bin_info$aln_number[2:length(medians_zoom_bin_info$aln_number)]
  mid_point <- (medians_zoom_bin_info$bin_start[2:length(medians_zoom_bin_info$aln_number)]+medians_zoom_bin_info$bin_end[2:length(medians_zoom_bin_info$aln_number)])/2
  medians_zoom_block_sum <- block_summary(medians_zoom_bin_info, bin_size = bin_size, show_only_coords = FALSE, reference_name = reference_assembly, query_name = query_assembly)
  if (print_tables) { 
    print(paste0("BINS AT ", bin_size/1000000, "-MBP BIN SIZE"))
    print(medians_zoom_bin_info)
    print(paste0("BLOCK SUMMARY AT ", bin_size/1000000, "-MBP BIN SIZE"))
    print(medians_zoom_block_sum)
  }
  medians_zoom_block_coords <- block_summary(medians_zoom_bin_info, bin_size = bin_size, show_only_coords = TRUE, reference_name = reference_assembly, query_name = query_assembly)
  ifelse(prediction_text, haploblock_prediction_text <- paste0("HAPLOBLOCK PREDICTIONS AT ", bin_size/1e06, "-MBP BIN SIZE"), haploblock_prediction_text <- "")
  ifelse(aln_text, aln_text <- paste0(aln_number_per_bin, " aln"), aln_text <- "")
  coords = data.frame()
  if(nrow(medians_zoom_block_coords) == 0){
    coords <- data.frame(start = zoom_start , end = zoom_end)
    fill_predictions <- "white"
    haploblock_prediction_text <- ""
  } else {
    coords <- block_summary(medians_zoom_bin_info, bin_size = bin_size, show_only_coords = TRUE, reference_name = reference_assembly, query_name = query_assembly)
  }
  graph <- plot_aln_pid_and_length(aln_subset, xmin = zoom_start, ymin = ymin, xmax = zoom_end, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = x_label_gap, dot_size = dot_size)
  if (is.data.frame(highlighted_target) == FALSE){
    graph_common <- graph + geom_vline(xintercept = seq(zoom_start, zoom_end, by = bin_size)) + geom_rect(data = coords, inherit.aes = FALSE, aes(xmin = coords[[1]], xmax = coords[[2]], ymin = ymin+1/4*(100-ymin), ymax = 100), color = "transparent", fill = fill_predictions, alpha = 0.3) +  annotate ("text", x = coords[[1]]+(coords[[2]]-coords[[1]])/2, y = (ymin+1/4*(100-ymin)+2/12*(100-ymin)), label = haploblock_prediction_text, size = 8, fontface = "bold", colour = color_prediction_text) +  annotate("text", x = mid_point, y = (ymin+1/4*(100-ymin)+1/12*(100-ymin)), size = 5, label = aln_text)
    print(graph_common)
  } else {
    graph_onlyiftarget <- graph + geom_vline(xintercept = seq(zoom_start, zoom_end, by = bin_size)) + geom_rect(data = coords, inherit.aes = FALSE, aes(xmin = coords[[1]], xmax = coords[[2]], ymin = ymin+1/4*(100-ymin), ymax = 100), color = "transparent", fill = fill_predictions, alpha = 0.3) +  annotate ("text", x = coords[[1]]+(coords[[2]]-coords[[1]])/2, y = (ymin+1/4*(100-ymin)+2/12*(100-ymin)), label = haploblock_prediction_text, size = 8, fontface = "bold", colour = color_prediction_text) +  annotate("text", x = mid_point, y = (ymin+1/4*(100-ymin)+1/12*(100-ymin)), size = 5, label = aln_text) + geom_rect(data = highlighted_target, inherit.aes = FALSE, aes(xmin = highlighted_target[[1]], xmax = highlighted_target[[2]], ymin = ymin, ymax = (ymin+1/4*(100-ymin))), color = "transparent", fill = fill_target, alpha = 0.3) + annotate ("text", x = highlighted_target[[1]]+(highlighted_target[[2]] - highlighted_target[[1]])/2, y = (ymin+1/6*(100-ymin)), label = target_text, size = 6, fontface = "bold", colour = color_target_text) + annotate ("text", x = highlighted_target[[1]]+(highlighted_target[[2]] - highlighted_target[[1]])/2, y = (ymin+1/12*(100-ymin)), label = paste0(nrow(aln_subset[(aln_subset$r_mid>=highlighted_target[[1]])&(aln_subset$r_mid<=highlighted_target[[2]]),]), " aln"), size = 5) 
    print(graph_onlyiftarget)
  }
}


plot_bins_and_selected_region <- function(aln_subset = data.frame(), bin_size = bin_size, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "target region", fill_target = "orange", color_target_text = "black", fill_predictions = "green", color_prediction_text = "black",  ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000, print_tables = TRUE, prediction_text = TRUE){
  medians_zoom <- plot_line_bin_median(data = aln_subset, bin_size = bin_size, bin_start = zoom_start, bin_end = zoom_end, ymin = ymin, cut_off = cut_off, reference_name = reference_assembly, query_name = query_assembly, x_label_gap = x_label_gap)
  medians_zoom_bin_info <- assign_blocks_mummer(medians_zoom[["data"]], original_file = aln_subset)
  aln_number_per_bin <- medians_zoom_bin_info$aln_number[2:length(medians_zoom_bin_info$aln_number)]
  mid_point <- (medians_zoom_bin_info$bin_start[2:length(medians_zoom_bin_info$aln_number)]+medians_zoom_bin_info$bin_end[2:length(medians_zoom_bin_info$aln_number)])/2
  medians_zoom_block_sum <- block_summary(medians_zoom_bin_info, bin_size = bin_size, show_only_coords = FALSE, reference_name = reference_assembly, query_name = query_assembly)
  if (print_tables) { 
    print(paste0("BINS AT ", bin_size/1000000, "-MBP BIN SIZE"))
    print(medians_zoom_bin_info)
    print(paste0("BLOCK SUMMARY AT ", bin_size/1000000, "-MBP BIN SIZE"))
    print(medians_zoom_block_sum)
  }
  medians_zoom_block_coords <- block_summary(medians_zoom_bin_info, bin_size = bin_size, show_only_coords = TRUE, reference_name = reference_assembly, query_name = query_assembly)
  ifelse(prediction_text == TRUE, haploblock_prediction_text <- paste0("HAPLOBLOCK PREDICTIONS AT ", bin_size/1e06, "-MBP BIN SIZE"), haploblock_prediction_text <- "")
  coords = data.frame()
  if(nrow(medians_zoom_block_coords) == 0){
    coords <- data.frame(start = zoom_start , end = zoom_end)
    fill_predictions <- "lightgrey"
    haploblock_prediction_text <- ""
  } else {
    coords <- block_summary(medians_zoom_bin_info, bin_size = bin_size, show_only_coords = TRUE, reference_name = reference_assembly, query_name = query_assembly)
  }
  graph <- plot_line_bin_median(data = aln_subset, bin_start = zoom_start, bin_end = zoom_end, bin_size = bin_size, cut_off = cut_off, ymin = ymin, ymax = 100, reference_name = reference_assembly, query_name = query_assembly, x_label_gap = x_label_gap)
  if (is.data.frame(highlighted_target) == FALSE){
    graph_common <- graph + geom_vline(xintercept = seq(zoom_start, zoom_end, by = bin_size)) + geom_rect(data = coords, inherit.aes = FALSE, aes(xmin = coords[[1]], xmax = coords[[2]], ymin = ymin+1/4*(100-ymin), ymax = 100), color = "transparent", fill = fill_predictions, alpha = 0.3) +  annotate ("text", x = coords[[1]]+(coords[[2]]-coords[[1]])/2, y = (ymin+1/4*(100-ymin)+2/12*(100-ymin)), label = haploblock_prediction_text, size = 8, fontface = "bold", colour = color_prediction_text) 
    print(graph_common)
  } else {
    graph_onlyiftarget <- graph + geom_vline(xintercept = seq(zoom_start, zoom_end, by = bin_size)) + geom_rect(data = coords, inherit.aes = FALSE, aes(xmin = coords[[1]], xmax = coords[[2]], ymin = ymin+1/4*(100-ymin), ymax = 100), color = "transparent", fill = fill_predictions, alpha = 0.3) +  annotate ("text", x = coords[[1]]+(coords[[2]]-coords[[1]])/2, y = (ymin+1/4*(100-ymin)+2/12*(100-ymin)), label = haploblock_prediction_text, size = 8, fontface = "bold", colour = color_prediction_text) + geom_rect(data = highlighted_target, inherit.aes = FALSE, aes(xmin = highlighted_target[[1]], xmax = highlighted_target[[2]], ymin = ymin, ymax = (ymin+1/4*(100-ymin))), color = "transparent", fill = fill_target, alpha = 0.3) + annotate ("text", x = highlighted_target[[1]]+(highlighted_target[[2]] - highlighted_target[[1]])/2, y = (ymin+1/6*(100-ymin)), label = target_text, size = 6, fontface = "bold", colour = color_target_text) + annotate ("text", x = highlighted_target[[1]]+(highlighted_target[[2]] - highlighted_target[[1]])/2, y = (ymin+1/12*(100-ymin)), label = paste0(nrow(aln_subset[(aln_subset$r_mid>=highlighted_target[[1]])&(aln_subset$r_mid<=highlighted_target[[2]]),]), " aln"), size = 5) 
    print(graph_onlyiftarget)
  }
}

##########  2. IDENTIFY GENES WITHIN THE JUST-MAPPED HAPLOBLOCK ##########  

cs_id_filler <- function(data, library, chr, ref.var){
  ref.var.code <- c("ARI", "CHI", "JAG", "JUL", "LAC", "LDM", "MAC", "NOR", "TSP", "STA", "SYM")
  ref.var.names <- c("arinalrfor", "chinese_spring", "jagger", "julius", "lancer", "landmark", "mace", "norin61", "spelta", "stanley", "mattis")
  names(ref.var.code) <- ref.var.names
  library <- library[grepl(paste("^TraesCS", chr, sep = ""), cs_id$chinese.spring.source) & grepl(paste("^Traes", ref.var.code[ref.var], chr, sep = ""), 
                                                                                                  cs_id$X10wheat.project.geneID),]
  data$cs_id <- NA
  for (i in 1:nrow(data)) {
    query <- data$var_id[i]
    replacement <- library$chinese.spring.source[grepl(query, library$X10wheat.project.geneID)]
    ifelse(length(replacement)==1, data$cs_id[i] <- replacement, data$cs_id[i] <- paste("source is not in chr", chr, "  ", sep = ""))
  }
  data$cs_id
  data$cs_id <- gsub(".{2}$", "", data$cs_id)
  return(data)
}

##########  3. PROVE IF THE HAPLOBLOCK WAS CALLED BY GENE-BASED BLAST PAIRWISE ALIGNMENTS ##########  

read_pairwise_position <- function(blast_path_gz, gtf, outfile, write_table = TRUE){
  all_comp <- read.table(gzfile(blast_path_gz), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  head(all_comp)
  
  #get the refseq position of the genes
  gtf$transcript_id <- str_split_fixed(gtf$V9, ";", 2)[,1]
  gtf$transcript_id <- gsub("transcript_id ", "", gtf$transcript_id)
  
  gene_only <- gtf[gtf$V3 == "gene",]
  
  gene_positions <- gene_only[,c("transcript_id", "V1", "V4",  "V5", "V7")]
  colnames(gene_positions) <- c("transcript", "chr", "start", "end", "strand")
  gene_positions$transcript <- gsub("ID=", "", gene_positions$transcript)
  
  #now add positions to the pairwise comparison file
  
  all_comp_positions <- merge(all_comp, gene_positions, all.x = TRUE, all.y = FALSE)
  head(all_comp_positions)
  
  if(write_table == TRUE){
    write.table(all_comp_positions, file = outfile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  return(all_comp_positions)
}

edited_calculate_pid_windows <- function(aln_data, window_size = 20, drop_no = 2){
  aln_data <- aln_data[order(aln_data$cs_start),]
  aln <- unique(aln_data$aln_type)
  mean_pidents <- data.frame(start = numeric(), 
                             end = numeric(), 
                             cs_start = numeric(), 
                             cs_end = numeric(), 
                             var_start = numeric(), 
                             var_end = numeric(),
                             pident_mean = numeric(), 
                             aln_type = character(),
                             start_cs_transcript = character(),
                             end_cs_transcript = character(),
                             start_var_transcript = character(),
                             end_var_transcript = character()) 
  
  for(j in seq(1, nrow(aln_data)-window_size)){
    window_start = j
    cs_start = aln_data[window_start, "cs_start"]
    var_start = aln_data[window_start, "var_start"]
    start_cs_transcript = aln_data[window_start, "cs_transcript"]
    start_var_transcript = aln_data[window_start, "var_transcript"]
    window_end = (window_start + window_size) - 1
    cs_end = aln_data[window_end, "cs_end"]
    var_end = aln_data[window_end, "var_end"]
    end_cs_transcript = aln_data[window_end, "cs_transcript"]
    end_var_transcript = aln_data[window_end, "var_transcript"]
    pidents <- aln_data[c(window_start:window_end), "pident"]
    pidents <- pidents[order(pidents, decreasing = TRUE)]
    pidents_sub <- pidents[c(1:(length(pidents)-drop_no))]
    pident_mean <- mean(pidents_sub)
    aln_to_add = data.frame(start = window_start, 
                            end = window_end, 
                            cs_start = cs_start, 
                            cs_end = cs_end, 
                            var_start = var_start,
                            var_end = var_end,
                            pident_mean = pident_mean, 
                            aln_type = aln, 
                            start_cs_transcript = start_cs_transcript,
                            end_cs_transcript = end_cs_transcript,
                            start_var_transcript = start_var_transcript,
                            end_var_transcript = end_var_transcript)
    mean_pidents <- rbind(mean_pidents, aln_to_add)
  }
  return(mean_pidents)
}

assign_blocks <-function(mean_pidents){
  mean_pidents_copy <- mean_pidents
  mean_pidents_copy$block_no <- NA
  block_no = 1
  new_block <- "no"
  
  for (i in seq(1, nrow(mean_pidents_copy))){
    #print(i)
    if(mean_pidents_copy[i, "pident_mean"] < 100){
      print("less than 100%")
      mean_pidents_copy[i, "block_no"] <- NA
    } else if (mean_pidents_copy[i, "pident_mean"] == 100){
      if(new_block == "no"){
        mean_pidents_copy[i, "block_no"] <- block_no
      } else if (new_block == "yes"){
        if(mean_pidents_copy[i, "start_position"] < end_prev_block){
          block_no <- block_no-1
          mean_pidents_copy[i, "block_no"] <- block_no
        } else {
          mean_pidents_copy[i, "block_no"] <- block_no
        }
      }
      if(i == nrow(mean_pidents_copy)){
        print("end")
      } else if (mean_pidents_copy[i+1, "pident_mean"] == 100){
        block_no <- block_no
        new_block <- "no"
        print("same block")
      } else if (mean_pidents_copy[i+1, "pident_mean"] < 100){
        block_no <- block_no+1
        new_block <- "yes"
        print("new_block")
        end_prev_block <- mean_pidents_copy[i, "end_position"]
      }
    }
  }
  return(mean_pidents_copy)
}
