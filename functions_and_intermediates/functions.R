# -------------------------------------------------------------------------
# File: functions.R
# Pierce et al 2025
# Purpose: Functions for IntronArchitecture
# -------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
library(data.table)
library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(Hmisc)
library(reshape2)
library(randomForest)
library(Biostrings)
library(GenomicRanges)
library(data.table)
library(RColorBrewer)
library(Rsamtools)
library(writexl)
library(cowplot)

# Functions ---------------------------------------------------------------
gffGenes <- function(input_gff) {

  # Read the GFF file
  gff_data <- fread(input_gff, sep = "\t", header = FALSE, data.table = FALSE)

  # Filter rows where the third column (type) is "mRNA"
  allGenes <- gff_data[gff_data[, 3] == "mRNA", ]

  # Remove unnecessary columns (V2, V3, V6, V8)
  allGenes <- allGenes[, -c(2, 3, 6, 8)]

  # Rename columns
  colnames(allGenes) <- c("chr", "start", "stop", "strand", "gene")

  # Clean and filter data
  allGenes$chr <- gsub("Chr", "", allGenes$chr)  # Remove "Chr" prefix
  allGenes$gene <- gsub("ID=", "", allGenes$gene)  # Remove "ID=" prefix
  allGenes$gene <- gsub(";.*", "", allGenes$gene)  # Remove trailing metadata

  allGenes$geneLength <- allGenes$stop - allGenes$start + 1  # Calculate gene length

  # Return the processed table
  return(allGenes)
}

gffExons <- function(input_gff) {

  # Read the GFF file
  gff_data <- fread(input_gff, sep = "\t", header = FALSE, data.table = FALSE)

  # Filter rows where the third column (type) is "exon"
  allExons <- gff_data[gff_data[, 3] == "exon", ]

  # Remove unnecessary columns (V2, V3, V6, V8)
  allExons <- allExons[, -c(2, 3, 6, 8)]

  # Rename columns
  colnames(allExons) <- c("chr", "start", "stop", "strand", "gene")

  # Clean and filter data
  allExons$chr <- gsub("Chr", "", allExons$chr)  # Remove "Chr" prefix
  allExons$gene <- gsub("Parent=", "", allExons$gene)  # Remove "ID=" prefix
  allExons$gene <- gsub("ID=", "", allExons$gene)  # Remove "ID=" prefix
  allExons$gene <- gsub(";.*", "", allExons$gene)  # Remove trailing metadata
  allExons$gene <- gsub(".exon.*", "", allExons$gene)  # Remove trailing metadata

  # Extract exon feature and calculate exon length
  allExons$feature <- "Exon"
  allExons$exonLength <- allExons$stop - allExons$start + 1

  # Assign exon numbers within each gene, considering strand
  allExons <- allExons %>%
    group_by(gene) %>%
    arrange(gene, if_else(strand == "+", start, -start), .by_group = TRUE) %>%  # Sort by strand
    mutate(ExonNum = row_number()) %>%  # Add sequential exon number
    ungroup()

  return(allExons)
}

annotateIntrons <- function(input_gff) {

  exon_data<-gffExons(input_gff)
  # Ensure input is a data.table
  exon_data <- as.data.table(exon_data)

  # Process the data to find intron coordinates
  intron_data <- exon_data %>%
    arrange(gene, start) %>%
    group_by(gene) %>%
    mutate(
      next_exon_start = lead(start),
      intron_start = stop + 1,
      intron_stop = next_exon_start - 1
    ) %>%
    filter(!is.na(intron_stop))  # Remove entries where no subsequent exon exists

  # Convert to data.table for further processing
  intron_data <- as.data.table(intron_data)

  # Select relevant columns for intron data
  intron_data <- intron_data[, .(gene, chr, intron_start, intron_stop, strand)]
  setnames(intron_data, c("gene", "chr", "intron_start", "intron_stop", "strand"), c("gene", "chr", "start", "stop", "strand"))

  # Calculate intron lengths
  intron_data[, intronLength := stop - start + 1]

  # Assign intron numbers based on order
  intron_data <- intron_data %>%
    group_by(gene) %>%
    arrange(gene, if_else(strand == "+", start, -start), .by_group = TRUE) %>%
    mutate(IntronNum = row_number()) %>%
    ungroup()

  intron_data$feature<-"Intron"

  # Return the annotated intron data
  return(as.data.table(intron_data))
}

intronFeatures <- function(input_gff) {

  allIntrons<-annotateIntrons(input_gff)
  allExons<-gffExons(input_gff)
  allGenes<-gffGenes(input_gff)

  allIntrons <- allIntrons %>% group_by(gene) %>%
    mutate(totalIntronLength=sum(intronLength))
  allIntrons<-as.data.table(allIntrons)

  firstIntrons<-allIntrons[IntronNum==1] %>% select(gene, intronLength, totalIntronLength)

  allExons<-left_join(allExons, firstIntrons, by="gene")
  allExons<-allExons %>% group_by(gene) %>%
    mutate(totalexonNum=n(), totalintronNum=n()-1, totalExonLength=sum(exonLength))
  allExons<-as.data.table(allExons)

  firstExons<-allExons[ExonNum==1]
  firstExons$firstintronpos<-ifelse(firstExons$totalintronNum == 0, NA, firstExons$exonLength + 1)
  names(firstExons)[names(firstExons) == "intronLength"] <- "firstIntronLength"
  names(firstExons)[names(firstExons) == "exonLength"] <- "firstExonLength"

  allGenes<-as.data.table(allGenes)
  firstExons <- left_join(firstExons, allGenes[, .(gene, geneLength)], by = "gene")

  firstExons <- firstExons[, !c("feature", "start","stop"), with = FALSE]

  return(list(allGenes=allGenes, allExons=allExons,allIntrons=allIntrons, IntronFeatures=firstExons))
}

isRepresentative <- function(data_table, representative_genes_file) {

  # Read the representative genes file
  representative_genes <- fread(representative_genes_file, header = FALSE, col.names = "gene")

  # Ensure data_table is a data.table
  data_table <- as.data.table(data_table)

  # Filter data_table to keep only rows with genes in the representative list
  filtered_data <- data_table[gene %in% representative_genes$gene]

  return(filtered_data)
}

isNuclear <- function(data_table) {
  # Keep only rows with numeric chromosome values
  data_table <- data_table[grepl("^[0-9]+$", data_table$chr), ]

  return(data_table)
}

firstIntron_density <- function(intron_features) {
  # Ensure input tables are data.tables
  merged_data <- as.data.table(intron_features)

  # Calculate relative intron start position
  merged_data[, rel_start := (firstExonLength + 1) / geneLength]

  # Filter out invalid or missing values
  merged_data <- merged_data[
    !is.na(rel_start) & !is.na(firstintronpos)
  ]

  # Create a density plot
  ggplot(merged_data[totalintronNum!=0], aes(x = rel_start)) +
    geom_density(fill = "#6a51a3", alpha = 0.8) +
    labs(
      x = "Relative Position in Gene Body",
      y = "Density",
      title = "Density of First Intron Positions Across Gene Bodies"
    ) +
    theme_classic(base_size = 6)
}

firstIntron_heatmap <- function(intron_features) {
  # Ensure input is a data.table
  merged_data <- as.data.table(intron_features)

  # Calculate relative intron start position
  merged_data[, rel_start := (firstExonLength + 1) / geneLength]

  # Filter out invalid or missing values
  merged_data <- merged_data[
    !is.na(rel_start) & !is.na(firstintronpos)
  ]

  # Create a 2D density heatmap
  ggplot(merged_data[totalintronNum!=0], aes(x = rel_start, y = 1)) +
    stat_density_2d(
      aes(fill = after_stat(density)),  # Updated to use after_stat
      geom = "tile",
      contour = FALSE,
      h = c(0.05, 1)
    ) +
    scale_fill_gradient(low = "#dadaeb", high = "#6a51a3") +
    labs(
      x = "Relative Position in Gene Body",
      y = "",
      title = "Density Heatmap of First Intron Positions",
      fill = "Density"
    ) +
    theme_minimal(base_size = 6) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
}

read_bw <- function(file){
  temp<-as.data.table(import.bw(file))
  colnames(temp)<-c("chr","start","stop", "width", "strand", "score")
  temp<-temp[,c("chr", "start", "stop", "score")]
  temp$chr<-gsub("chr", "", temp$chr)
  temp$chr<-as.character(temp$chr)
  return(temp)
}

chrom10_make_windows <- function(chromosome, loc_start, loc_end, windowsize=10){
  window_table<-data.table(start=seq(loc_start,loc_end-10,by=10))
  window_table$stop<-seq(loc_start+10,loc_end, by=10)
  window_table<-cbind(chr = chromosome, window_table)

  window_table$window_id<-1:nrow(window_table)
  setkey(window_table, chr, start, stop)
}

overlaps_bw <- function(temp, data_table){
  overlaps<-foverlaps(temp, data_table, by.x = c("chr", "start", "stop"),
                      type="any")
  overlaps[, len := pmax(0, pmin(stop, i.stop) - pmax(start, i.start) + 1)]
  overlaps[, totalDepth := score*len]
  result <- overlaps[, .(chr=first(chr), start=first(start), stop=first(stop),
                         sum_overlaps = sum(totalDepth)), by = .(window_id)]
  result<-na.omit(result)

  return(result)
}

overlaps_files_list <- function(data_table, files_list){
  dt<-data_table
  for (i in files_list) {
    message(i)
    columnname<-gsub("\\.bw","", i)
    columnname<-gsub(".*/", "", columnname)
    temp<-read_bw(i)
    overlaps<-overlaps_bw(temp, data_table)
    dt$tempname<-overlaps$sum_overlaps[match(dt$window_id, overlaps$window_id)]
    names(dt)[names(dt) == "tempname"] <- columnname
  }

  return(dt)
}

epigenomeOverlaps <- function(organism_folder_path, filter_chr = NULL) {
  # List all files in the organism's folder
  bw_files <- list.files(organism_folder_path, pattern = "\\.bw$", full.names = TRUE)

  # Load BW file and convert to data table
  bw_data <- read_bw(bw_files[1])

  # Filter chromosomes if filter_chr is specified
  if (!is.null(filter_chr)) {
    bw_data <- bw_data[grep(paste0("^", filter_chr), chr, ignore.case = TRUE)]
  }

  # Identify unique chromosomes and their maximum lengths
  chromosomes_info <- bw_data[, .(max_length = max(stop)), by = chr]

  # Initialize a list to store results
  scores_list <- list()

  # Process each chromosome
  for (i in seq_len(nrow(chromosomes_info))) {
    chr <- chromosomes_info$chr[i]
    max_length <- chromosomes_info$max_length[i]

    # Inform which chromosome is being processed
    message("Processing chromosome: ", chr)

    # Generate windows and calculate scores
    chrom_data <- chrom10_make_windows(chr, 0, max_length)
    chrom_data$chr <- as.character(chrom_data$chr)
    setkey(chrom_data, chr, start, stop)

    scores <- overlaps_files_list(chrom_data, bw_files)
    scores_list[[i]] <- scores
  }

  # Combine results for all chromosomes
  all_scores <- rbindlist(scores_list)

  return(all_scores)
}

featureOverlaps <- function(organism_folder_path, windows) {
  # List all files in the organism's folder
  bw_files <- list.files(organism_folder_path, pattern = "\\.bw$", full.names = TRUE)

  setkey(windows, chr, start, stop)

  scores <- overlaps_files_list(windows, bw_files)

  return(scores)
}

IntronExonScores <- function(exon_coords, intron_coords, features, pattern = "") {
  # Function to process overlaps and aggregate scores
  setkey(exon_coords, chr, start, stop)
  setkey(intron_coords, chr, start, stop)
  setkey(features, chr, start, stop)
  process_overlap <- function(coords, features, id_col, length_col, feature_label) {
    # Find overlaps
    overlap_results <- foverlaps(coords, features, by.x = c("chr", "start", "stop"), type = "any")
    overlap_results <- as.data.table(overlap_results)

    # Aggregate sums based on mean
    aggregated_sums <- overlap_results[, lapply(.SD, mean, na.rm = TRUE),
                                       by = .(gene, get(id_col)),
                                       .SDcols = patterns(pattern)]

    # Rename dynamically evaluated column for consistency
    setnames(aggregated_sums, "get", id_col)

    # Merge aggregated results back with original coordinates
    coords <- merge(coords, aggregated_sums, by = c("gene", id_col), all.x = TRUE)

    # Add feature label and rename columns
    coords$Feature <- feature_label
    setnames(coords, id_col, "Num")
    setnames(coords, length_col, "Length")

    return(coords)
  }

  # Process exon overlaps
  exon_coords <- process_overlap(exon_coords, features, id_col = "ExonNum", length_col = "exonLength", feature_label = "Exon")

  # Process intron overlaps
  intron_coords <- process_overlap(intron_coords, features, id_col = "IntronNum", length_col = "intronLength", feature_label = "Intron")

  # Combine exon and intron results
  final_df <- rbind(exon_coords, intron_coords, fill = TRUE)

  return(final_df)
}

geneScores <- function(gene_coords, features, pattern = "") {
  # Function to process overlaps and aggregate scores
  setkey(gene_coords, chr, start, stop)
  setkey(features, chr, start, stop)
  # Find overlaps
  overlap_results <- foverlaps(gene_coords, features, by.x = c("chr", "start", "stop"), type = "any")
  overlap_results <- as.data.table(overlap_results)

  # Aggregate sums based on mean
  aggregated_sums <- overlap_results[, lapply(.SD, mean, na.rm = TRUE),
                                     by = .(gene),
                                     .SDcols = patterns(pattern)]

  # Merge aggregated results back with original coordinates
  gene_coords <- merge(gene_coords, aggregated_sums, by = "gene", all.x = TRUE)

  return(gene_coords)
}

summarizeMarksbyGroup <- function(db, epigmarks, groupingVar) {
  summary <- db %>%
    pivot_longer(cols = all_of(epigmarks), names_to = "Mark", values_to = "Value") %>%
    group_by(across(all_of(groupingVar)), Mark, Feature) %>%
    summarise(`Number of Genes` = n(),
              Mean = mean(Value, na.rm = TRUE),
              SE = sd(Value, na.rm = TRUE) / sqrt(n()), .groups = 'drop')

  summary <- as.data.table(summary)
  summary$Feature <- as.factor(summary$Feature)

  return(summary)
}
plotSplitData <- function(split_data) {
  plotlist <- list()  # Initialize an empty list

  for (mark in names(split_data)) {
    data <- split_data[[mark]]  # Extract data for the current mark

    # Create the plot
    plot <- ggplot(data, aes(x = featureNum, y = Mean, fill = Feature)) +
      geom_point(shape = 21, col = "black") +
      theme_classic(base_size = 6) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      geom_linerange(
        aes(ymin = Mean - SE, ymax = Mean + SE),
        color = "black",
        alpha = 0.6
      ) +
      scale_fill_manual(values = c("#74add1", "#fdae80")) +
      labs(y = mark)

    # Add the plot to the list using double brackets
    plotlist[[mark]] <- plot
  }

  return(plotlist)  # Return the list of plots
}

plotSplitDataV2 <- function(split_data) {
  plotlist <- list()  # Initialize an empty list

  for (mark in names(split_data)) {
    data <- split_data[[mark]]  # Extract data for the current mark

    # Convert featureNum to numeric (if it's a factor)
    data$featureNum <- as.numeric(as.character(data$featureNum))

    # Identify positions where Feature contains "Intron"
    intron_positions <- data$featureNum[grepl("Intron", data$Feature)]

    # Create a data frame to define shaded regions
    shading_df <- data.frame(
      xmin = intron_positions - 0.5,
      xmax = intron_positions + 0.5,
      ymin = -Inf,
      ymax = Inf
    )

    # Create the plot
    plot <- ggplot(data, aes(x = featureNum, y = Mean, fill = Feature)) +
      # Shading for intron regions
      geom_rect(data = shading_df,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = "grey90", alpha = 0.5) +
      geom_point(shape = 21, col = "black") +
      geom_linerange(aes(ymin = Mean - SE, ymax = Mean + SE),
                     color = "black", alpha = 0.6) +
      theme_classic(base_size = 6) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      scale_fill_manual(values = c("#74add1", "#fdae80")) +
      labs(y = mark)

    # Add the plot to the list
    plotlist[[mark]] <- plot
  }

  return(plotlist)  # Return the list of plots
}

calculate_enrichment <- function(data, mark_cols, input_col) {
  # Hardcoded column names for start and stop
  start_col <- "start"
  stop_col <- "stop"

  # Calculate total depth for input once
  input_depth <- sum(data[[input_col]] * (data[[stop_col]] - data[[start_col]]), na.rm = TRUE)

  # Loop through each mark column
  for (mark_col in mark_cols) {
    # Calculate total depth for the current mark
    mark_depth <- sum(data[[mark_col]] * (data[[stop_col]] - data[[start_col]]), na.rm = TRUE)

    # Add enrichment values as a new column
    enrichment_col <- paste0(mark_col, "_enrichment")
    data[[enrichment_col]] <- log2((1 + data[[mark_col]]) / mark_depth) -
      log2((1 + data[[input_col]]) / input_depth)
  }

  # Return the updated data table
  return(data)
}

library(data.table)

# Define the function
file_correlations <- function(data_table, excluded_cols = c("chr", "start", "stop", "window_id"), output_dir = "./") {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Remove excluded columns
  data <- data_table[, !..excluded_cols, with = FALSE]

  # Group columns by prefix
  get_prefix <- function(colname) {
    strsplit(colname, "-")[[1]][1]
  }
  column_groups <- split(names(data), sapply(names(data), get_prefix))

  # Initialize lists to store correlation results
  spearman_results <- list()
  pearson_results <- list()

  # Calculate pairwise correlations for each group
  for (group in names(column_groups)) {
    group_cols <- column_groups[[group]]

    # Subset the data for the group
    group_data <- data[, ..group_cols]

    # Calculate Spearman and Pearson correlations
    spearman_corr <- cor(group_data, method = "spearman", use = "pairwise.complete.obs")
    pearson_corr <- cor(group_data, method = "pearson", use = "pairwise.complete.obs")

    # Store the results
    spearman_results[[group]] <- spearman_corr
    pearson_results[[group]] <- pearson_corr

    # Export to CSV
    fwrite(as.data.table(spearman_corr, keep.rownames = "Variable"),
           file = file.path(output_dir, paste0("spearman_", group, "_correlations.csv")))
    fwrite(as.data.table(pearson_corr, keep.rownames = "Variable"),
           file = file.path(output_dir, paste0("pearson_", group, "_correlations.csv")))
  }

  print("Correlation files saved for each group.")
  return(list(spearman = spearman_results, pearson = pearson_results))
}

library(data.table)


normalizedSum <- function(data_table, excluded_cols = c("chr", "start", "stop", "window_id")) {
  # Group columns by prefix
  get_prefix <- function(colname) {
    strsplit(colname, "-")[[1]][1]
  }

  # Retain excluded columns for the final output
  retained_data <- data_table[, ..excluded_cols, with = FALSE]

  # Exclude columns for normalization calculation
  data_table <- data_table[, !..excluded_cols, with = FALSE]
  column_groups <- split(names(data_table), sapply(names(data_table), get_prefix))

  # Initialize data.table to store results
  normalized_data <- data.table()

  # Calculate normalized sum for each group
  for (group in names(column_groups)) {
    group_cols <- column_groups[[group]]

    # Subset the data for the group
    group_data <- data_table[, ..group_cols]

    # Calculate sum and count non-NA values
    sum_col <- rowSums(group_data, na.rm = TRUE)
    col_count <- rowSums(!is.na(group_data))

    # Calculate normalized sum
    norm_sum <- ifelse(col_count > 0, sum_col / col_count, NA)

    # Add normalized sum to the data.table
    normalized_data <- cbind(normalized_data, setNames(data.table(norm_sum), paste0("NormalizedSum_", group)))
  }

  # Combine retained columns with normalized data
  final_data <- cbind(retained_data, normalized_data)

  print("Normalized sums calculated for each group.")
  return(final_data)
}

correlationDataTable <- function(data_table, var1_columns, var2_columns) {
  # Calculate Spearman correlation matrix
  cortable <- cor(data_table, use = "complete.obs", method = "spearman")

  # Melt the correlation table for heatmap visualization
  long_cortable <- reshape2::melt(cortable)
  long_cortable$Var1 <- as.character(long_cortable$Var1)
  long_cortable$Var2 <- as.character(long_cortable$Var2)

  # Filter Var1 columns
  filtertable <- dplyr::filter(long_cortable, Var1 %in% var1_columns)

  # Filter Var2 columns
  filtertable <- dplyr::filter(filtertable, Var2 %in% var2_columns)

  return(filtertable)
}

rankColumns <- function(data_table, columns_to_rank, suffix = "_ranked", ties_method = "average") {
  # Ensure all columns to rank exist in the data frame
  columns_to_rank <- intersect(columns_to_rank, colnames(data_table))

  # Apply ranking to specified columns and append ranked columns to the data frame
  for (col in columns_to_rank) {
    ranked_col_name <- paste0(col, suffix)
    data_table[[ranked_col_name]] <- rank(data_table[[col]], ties.method = ties_method)
  }

  # Return the modified data frame
  return(data_table)
}

combineMetrics <- function(mark, CV_list, j) {
  # Initialize an empty list to store the metrics
  metrics_list <- list()

  # Loop through the CV_list
  for (i in seq_along(CV_list)) {
    metrics_list[[i]] <- do.call(rbind, CV_list[[i]][, j])
  }

  # Combine the list of metrics into a single data frame
  metrics <- do.call(rbind, metrics_list)

  # Add the "mark" column
  metrics$mark <- mark

  return(metrics)
}

RFMarks <- function(data, marks, predictors, ntree = 150) {
  models <- list()

  for (mark in marks) {
    formula_string <- paste(mark, "~", paste(predictors, collapse = "+"))
    rf_formula <- as.formula(formula_string)

    model <- randomForest(rf_formula, data = data, ntree = ntree, importance = TRUE)
    models[[mark]] <- model
  }

  return(models)
}

perform_cv_once <- function(i, df, response_variable) {
  formula <- as.formula(paste(response_variable,
                              "~ totalintronNum + firstIntronLength + firstintronpos + totalIntronLength + geneLength + IntronPct + totalExonLength"))
  train1 <- df[rand != i]
  test1 <- df[rand == i]
  rf_train <- randomForest(formula,
                           data = train1, ntree = 150, importance = TRUE)
  featureimp<-as.data.table(importance(rf_train), keep.rownames = TRUE)
  featureimp$CV<-i

  test1$predicted <- predict(rf_train, test1)

  RMSE <- sqrt(mean((test1$predicted - test1[[response_variable]])^2))
  MAE <- mean(abs(test1$predicted - test1[[response_variable]]))
  R2 <- cor(test1$predicted, test1[[response_variable]])^2
  COR <- cor(test1$predicted, test1[[response_variable]])
  return(list(data.table(test1,i, RMSE, MAE, R2, COR), featureimp))

}

# Function to repeat the CV process n times and store results in a list
repeat_cv_n_times <- function(df, n, response_variable) {
  list_of_results <- vector("list", n)
  list_of_results2 <- vector("list", n)

  for(i in 1:n) {
    df$rand<-sample(rep(1:5, length.out = nrow(df)))

    for (j in 1:5) {
      list_of_results[[j]] <- perform_cv_once(j, df, response_variable)
    }
    # Combine the results into a single data frame for this iteration
    combined_results <- do.call(rbind, list_of_results)

    # Append the combined results to the list
    list_of_results2[[i]] <- combined_results
  }

  return(list_of_results2)
}

processCrossvalidation <- function(data, column_names, n_times) {
  results <- lapply(column_names, function(column_name) {
    # Perform cross-validation
    cv_result <- repeat_cv_n_times(data, n_times, column_name)

    # Extract metrics and features
    metrics <- combineMetrics(column_name, cv_result, 1)
    features <- combineMetrics(column_name, cv_result, 2)

    # Return a list of results for this column
    list(
      CV = cv_result,
      Metrics = metrics,
      Features = features
    )
  })

  # Assign column names to the results list for clarity
  names(results) <- column_names
  return(results)
}

introngroup_summaries <- function(data, marks,intronGroup=intronGroup) {
  summary <- data %>%
    pivot_longer(cols = all_of(marks), names_to = "Mark", values_to = "Value") %>%
    group_by(intronGroup, Mark) %>%
    summarise(`Number of Genes` = n(),
              Mean = mean(Value, na.rm = TRUE),
              SE = sd(Value, na.rm = TRUE) / sqrt(n()))

  return(summary)
}

plot_introngroup2<- function(summary, mark, genefeature){
  plot<-ggplot(all[all$Mark==mark,], aes(x=intronGroup, y=Mean))+
    geom_point(shape=21, size=1.8, fill="#cc4c02")+
    geom_linerange(aes(ymin = Mean - SE,
                       ymax = Mean + SE), color="black", alpha=0.6)+
    theme_classic(base_size = 6)+
    #scale_fill_viridis(option="mako", begin = 0.4,end=0.9)+
    #scale_color_viridis(option="mako", begin = 0.4,end=0.9)+
    labs(x=genefeature, y=mark)+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  return(plot)

}

# From Polymorphology2 github.com/greymonroe/polymorphology2
plot_grid2 <- function(plotlist, type = "cols", relsize = rep(1, length(plotlist))) {
  if (!type %in% c("cols", "rows")) {
    stop("Invalid 'type'. Must be 'cols' or 'rows'.")
  }

  if (length(plotlist) != length(relsize)) {
    stop("Length of 'relsize' must match the number of plots in 'plotlist'.")
  }

  # Convert ggplot objects to grobs
  grobs <- lapply(plotlist, ggplotGrob)

  # Adjust widths or heights based on the type
  if (type == "cols") {
    for (i in seq_along(grobs)) {
      grobs[[i]]$widths <- grid::unit.pmax(grobs[[i]]$widths, grid::unit(relsize[i], "null"))
    }
    combined_grob <- do.call(cbind, grobs)
  } else if (type == "rows") {
    for (i in seq_along(grobs)) {
      grobs[[i]]$heights <- grid::unit.pmax(grobs[[i]]$heights, grid::unit(relsize[i], "null"))
    }
    combined_grob <- do.call(rbind, grobs)
  }

  # Draw the combined grob
  #grid::grid.newpage()
  grid::grid.draw(combined_grob)
}

extractCytme <- function(fasta_index, allGenes, methylation_data) {
  # Create GRanges object for the input data table
  regions <- GRanges(
    seqnames = Rle(allGenes$chr),
    ranges = IRanges(start = allGenes$start, end = allGenes$stop),
    gene = allGenes$gene
  )

  # Extract sequences from the FASTA file
  extracted_sequences <- getSeq(fasta_index, regions)
  names(extracted_sequences) <- mcols(regions)$gene

  # Ensure names align properly
  if (!all(names(extracted_sequences) == mcols(regions)$gene)) {
    stop("Gene names do not match between sequences and regions.")
  }

  # Count cytosines in each sequence
  cytosine_counts <- sapply(extracted_sequences, function(seq) sum(letterFrequency(seq, letters = "C", as.prob = FALSE)))
  gene_cytosine_counts <- data.frame(
    gene = names(extracted_sequences),
    Cytosine_Count = cytosine_counts
  )

  # Join cytosine counts with gene metadata
  all <- inner_join(gene_cytosine_counts, allGenes, by = "gene")
  all <- as.data.table(all)
  setkey(all, chr, strand, start, stop)

  # Find overlaps with methylation data
  overlaps <- foverlaps(all, methylation_data, by.x = c("chr", "strand", "start", "stop"), type = "any")
  overlaps[, c("V5", "V6") := NULL]

  # Group and summarize overlaps
  meOverlaps <- overlaps %>%
    group_by(gene, Context) %>%
    reframe(
      gene = gene,
      start = i.start,
      stop = i.stop,
      strand = strand,
      geneLength = geneLength,
      CytCount = Cytosine_Count,
      meCount = sum(Value, na.rm = TRUE)
    )

  meOverlaps <- unique(meOverlaps)
  meOverlaps$meCount[is.na(meOverlaps$meCount)] <- 0
  meOverlaps$ratio <- meOverlaps$meCount / meOverlaps$CytCount

  # Pivot the data to a wide format
  wide <- meOverlaps %>%
    pivot_wider(
      id_cols = c(gene, start, stop, strand, geneLength, CytCount),
      names_from = Context,
      values_from = ratio
    )

  # Convert to data.table and handle missing values
  wide <- as.data.table(wide)
  wide[, `NA` := NULL]
  wide$CG[is.na(wide$CG)] <- 0
  wide$CHH[is.na(wide$CHH)] <- 0
  wide$CHG[is.na(wide$CHG)] <- 0

  return(wide)
}
intronFeatureScatterPlots <- function(data, column_name, epigmarks, g = 10) {
  library(ggplot2)

  # Step 1: Create groups based on the specified column
  group_column_name <- paste0(column_name, "Group")
  data$intronGroup <- cut2(data[[column_name]], g = g)

  # Step 2: Summarize data by groups
  all_summaries <- introngroup_summaries(data, epigmarks)

  # Step 3: Generate plots for each epigenetic mark
  plots_list <- lapply(epigmarks, function(mark) {
    plot_introngroup2(all_summaries, mark, column_name)
  })

  # Step 4: Adjust plots
  plots_list <- lapply(seq_along(plots_list), function(i) {
    if (i <= length(plots_list) - 4) {
      # Remove x-axis text and title for all but the last 4 plots
      plots_list[[i]] + theme(axis.text.x = element_blank(),
                              axis.title.x = element_blank())
    } else {
      # Keep x-axis for the last 4 plots
      plots_list[[i]]
    }
  })

  plot_grid(plotlist = plots_list, ncol = 4, align = "v", rel_heights = 1)
}

metricsViolinPlot <- function(data, y_axis, x_axis = "mark", base_size = 6) {
  library(ggplot2)

  ggplot(data, aes(x = reorder(!!sym(x_axis), !!sym(y_axis)), y = !!sym(y_axis))) +
    geom_violin() +
    geom_jitter(color="#0570b0", alpha=0.4,size = 0.35, shape=16) +
    #geom_boxplot(width = 0.2, color = "grey40", alpha = 0.4) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      axis.title.x = element_blank()
    )
}
summarizeData <- function(db, epigmarks, groupingVar) {
  summary <- db %>%
    pivot_longer(cols = all_of(epigmarks), names_to = "Mark", values_to = "Value") %>%
    group_by(across(all_of(groupingVar)), Mark, n) %>%
    summarise(`Number of Genes` = n(),
              Mean = mean(Value, na.rm = TRUE),
              SE = sd(Value, na.rm = TRUE) / sqrt(n()), .groups = 'drop')

  summary <- as.data.table(summary)
  #summary$Feature <- as.factor(summary$Feature)

  return(summary)
}

featuresMetaplots <- function(data, group_by_column, grouping_column = "intronGp", marks, g = 10) {
  # Step 1: Create groups using cut2
  data[[grouping_column]] <- cut2(data[[group_by_column]], g = g)

  # Step 2: Summarize data (replace summarizeData with your actual function or summarize manually)
  summary <- summarizeData(data, marks, grouping_column)

  # Step 3: Create plots for each mark
  plot_list <- lapply(marks, function(mark) {
    ggplot(summary[summary$Mark == mark, ], aes(x = n, y = Mean, col = !!sym(grouping_column), group = !!sym(grouping_column))) +
      geom_line(linewidth = 0.5) +
      theme_classic(base_size = 6) +
      scale_color_manual(values = c("#d0d1e9", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#045a8d", "#08306B","#023851", "grey15", "black")) +
      #xlim(0, 300) +
      theme(
        legend.key.size = unit(0.3, "cm"),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(y = mark)+
      scale_x_continuous(limits = c(0, 300),
                         breaks = c(0, 100,150, 200, 300),
                         labels = c("-1500", "-500", "TSS", "500", "1500"))
  })

  names(plot_list) <- marks  # Name the plots by marks
  return(plot_list)
}
featuresGeneMetaplots <- function(data, group_by_column, grouping_column = "intronGp", marks, g = 10) {
  # Step 1: Create groups using cut2
  data[[grouping_column]] <- cut2(data[[group_by_column]], g = g)

  # Step 2: Summarize data (replace summarizeData with your actual function or summarize manually)
  summary <- summarizeData(data, marks, grouping_column)

  # Step 3: Create plots for each mark
  plot_list <- lapply(marks, function(mark) {
    ggplot(summary[summary$Mark == mark, ], aes(x = n, y = Mean, col = !!sym(grouping_column), group = !!sym(grouping_column))) +
      geom_line(linewidth = 0.5) +
      theme_classic(base_size = 6) +
      scale_color_manual(values = c("#d0d1e9", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#045a8d", "#08306B","#023851", "grey15", "black")) +
      #xlim(0, 600) +
      theme(
        legend.key.size = unit(0.3, "cm"),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(y = mark)+
      scale_x_continuous(limits = c(0, 600),
                         breaks = c(0, 200,300, 400, 600),
                         labels = c("-3000", "-1000", "TSS", "1000", "3000"))
  })

  names(plot_list) <- marks  # Name the plots by marks
  return(plot_list)
}
intron_marks_dups<-function(dups, gene_data, mark){

  dups$mark.1<-unlist(gene_data[,mark, with=F])[match(dups$dup1, gene_data$gene)]
  dups$mark.2<-unlist(gene_data[,mark, with=F])[match(dups$dup2, gene_data$gene)]
  dups$mark_diff<-dups$mark.1-dups$mark.2

  abs<-dups

  abs$mark_diff[abs$loss=="loss"]<- -abs$mark_diff[abs$loss=="loss"] #made "loss" ones negative because we dont know if it was loss or gain. But we have to subtract one column from the other

  sum<-abs[!is.na(mark_diff),.(mark_diff=mean(mark_diff),
                               mark_diff_se=2*sd(mark_diff)/sqrt(.N),
                               N=.N,
                               meanintrondiff=mean(abs(intron_diff)) ), by=.(intron=(cut(abs(intron_diff), breaks=c(-.1,.1, 1,3,Inf), labels=c("0", "+1", "+2,+3", "+>3"))))]
  sum$mark<-mark

  sum2<-abs[!is.na(mark_diff),.(mark_diff=mean(mark_diff),
                                mark_diff_se=2*sd(mark_diff)/sqrt(.N),
                                N=.N,
                                meanintrondiff=mean(abs(intron_diff)) ), by=.(intron=(cut(abs(intron_diff), breaks=c(-.1,.1,Inf), labels=c("0", "+"))))]
  sum2$mark<-mark

  ttest<-t.test(abs$mark_diff~abs(abs$intron_diff)>0)


  plot<-ggplot(sum, aes(x=intron, y=mark_diff,fill=intron))+
    geom_bar(stat="identity", width=0.8, col="black", linewidth=0.4)+
    scale_fill_manual(values=c("#f1eef6","#bdc9e1","#74a9cf","#0570b0"))+
    geom_errorbar(aes(ymin=mark_diff-mark_diff_se, ymax=mark_diff+mark_diff_se), width=0.2)+
    theme_classic(base_size = 6)+
    #geom_hline(yintercept = 0)+
    scale_y_continuous(name=mark)+
    labs(x="intron difference")+
    theme(legend.position = "none",
          axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          axis.ticks = element_blank())
  return(list(data=sum,data2=sum2, ttest=ttest, plot=plot))
}

intron_marks_transposed<-function(dups, gene_data, mark){

  dups$mark.1<-unlist(gene_data[,mark, with=F])[match(dups$Transposed, gene_data$gene)]
  dups$mark.2<-unlist(gene_data[,mark, with=F])[match(dups$Parental, gene_data$gene)]
  dups$mark_diff<-dups$mark.1-dups$mark.2

  abs<-dups

  sum<-abs[!is.na(mark_diff),.(mark_diff=mean(mark_diff),
                               mark_diff_se=2*sd(mark_diff)/sqrt(.N),
                               N=.N,
                               meanintrondiff=mean(intron_diff) ),
           by=.(intron=loss)]
  sum$mark<-mark
  sum$intron=factor(sum$intron, levels=c("loss", "equal","gain"))
  abs$intron=factor(abs$loss, levels=c("loss", "equal","gain"))
  ttest<-t.test(abs[intron!="equal"]$mark_diff~abs[intron!="equal"]$intron)


  plot<-ggplot(sum, aes(x=intron, y=mark_diff,fill=intron))+
    geom_bar(stat="identity", width=0.8, col="black", linewidth=0.4)+
    scale_fill_manual(values=c("#fee5d9","#fcaf99","#cc4c02"))+
    geom_errorbar(aes(ymin=mark_diff-mark_diff_se, ymax=mark_diff+mark_diff_se), width=0.2)+
    theme_classic(base_size = 6)+
    #geom_hline(yintercept = 0)+
    scale_y_continuous(name=mark)+
    labs(x="intron difference")+
    theme(legend.position = "none",
          axis.text.x=element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          axis.ticks = element_blank())
  return(list(data=sum, plot=plot, ttest=ttest))
}
extractCytmeIntron <- function(fasta_index, allGenes, methylation_data) {
  # Create GRanges object for the input data table
  regions <- GRanges(
    seqnames = Rle(allGenes$chr),
    ranges = IRanges(start = allGenes$start, end = allGenes$stop),
    gene = allGenes$gene,
    IntronNum = allGenes$IntronNum
  )

  # Extract sequences from the FASTA file
  extracted_sequences <- getSeq(fasta_index, regions)
  names(extracted_sequences) <- mcols(regions)$gene

  # Ensure names align properly
  if (!all(names(extracted_sequences) == mcols(regions)$gene)) {
    stop("Gene names do not match between sequences and regions.")
  }

  # Count cytosines in each sequence
  cytosine_counts <- sapply(extracted_sequences, function(seq) sum(letterFrequency(seq, letters = "C", as.prob = FALSE)))
  # gene_cytosine_counts <- data.frame(
  #   gene = names(extracted_sequences),
  #   Cytosine_Count = cytosine_counts
  #
  # )
  gene_cytosine_counts <- data.frame(
    gene = mcols(regions)$gene,
    IntronNum = mcols(regions)$IntronNum,
    Cytosine_Count = cytosine_counts
  )

  # Join cytosine counts with gene metadata
  all <- inner_join(gene_cytosine_counts, allGenes, by = c("gene", "IntronNum"))
  all <- as.data.table(all)
  setkey(all, chr, strand, start, stop)

  # Find overlaps with methylation data
  overlaps <- foverlaps(all, methylation_data, by.x = c("chr", "strand", "start", "stop"), type = "any")
  overlaps[, c("V5", "V6") := NULL]

  # Group and summarize overlaps
  meOverlaps <- overlaps %>%
    group_by(gene, Context, IntronNum) %>%
    reframe(
      gene = gene,
      start = i.start,
      stop = i.stop,
      strand = strand,
      IntronNum = IntronNum,
      Feature = feature,
      intronLength = intronLength,
      CytCount = Cytosine_Count,
      meCount = sum(Value, na.rm = TRUE)
    )

  meOverlaps <- unique(meOverlaps)
  meOverlaps$meCount[is.na(meOverlaps$meCount)] <- 0
  meOverlaps$ratio <- meOverlaps$meCount / meOverlaps$CytCount

  # Pivot the data to a wide format
  wide <- meOverlaps %>%
    pivot_wider(
      id_cols = c(gene, start, stop, strand, intronLength, CytCount),
      names_from = Context,
      values_from = ratio
    )

  # Convert to data.table and handle missing values
  wide <- as.data.table(wide)
  wide[, `NA` := NULL]
  wide$CG[is.na(wide$CG)] <- 0
  wide$CHH[is.na(wide$CHH)] <- 0
  wide$CHG[is.na(wide$CHG)] <- 0

  return(wide)
}
extractCytmeExon <- function(fasta_index, allGenes, methylation_data) {
  # Create GRanges object for the input data table
  regions <- GRanges(
    seqnames = Rle(allGenes$chr),
    ranges = IRanges(start = allGenes$start, end = allGenes$stop),
    gene = allGenes$gene,
    ExonNum = allGenes$ExonNum
  )

  # Extract sequences from the FASTA file
  extracted_sequences <- getSeq(fasta_index, regions)
  names(extracted_sequences) <- mcols(regions)$gene

  # Ensure names align properly
  if (!all(names(extracted_sequences) == mcols(regions)$gene)) {
    stop("Gene names do not match between sequences and regions.")
  }

  # Count cytosines in each sequence
  cytosine_counts <- sapply(extracted_sequences, function(seq) sum(letterFrequency(seq, letters = "C", as.prob = FALSE)))
  # gene_cytosine_counts <- data.frame(
  #   gene = names(extracted_sequences),
  #   Cytosine_Count = cytosine_counts
  #
  # )
  gene_cytosine_counts <- data.frame(
    gene = mcols(regions)$gene,
    ExonNum = mcols(regions)$ExonNum,
    Cytosine_Count = cytosine_counts
  )

  # Join cytosine counts with gene metadata
  all <- inner_join(gene_cytosine_counts, allGenes, by = c("gene", "ExonNum"))
  all <- as.data.table(all)
  setkey(all, chr, strand, start, stop)

  # Find overlaps with methylation data
  overlaps <- foverlaps(all, methylation_data, by.x = c("chr", "strand", "start", "stop"), type = "any")
  overlaps[, c("V5", "V6") := NULL]

  # Group and summarize overlaps
  meOverlaps <- overlaps %>%
    group_by(gene, Context, ExonNum) %>%
    reframe(
      gene = gene,
      start = i.start,
      stop = i.stop,
      strand = strand,
      ExonNum = ExonNum,
      Feature = feature,
      exonLength = exonLength,
      CytCount = Cytosine_Count,
      meCount = sum(Value, na.rm = TRUE)
    )

  meOverlaps <- unique(meOverlaps)
  meOverlaps$meCount[is.na(meOverlaps$meCount)] <- 0
  meOverlaps$ratio <- meOverlaps$meCount / meOverlaps$CytCount

  # Pivot the data to a wide format
  wide <- meOverlaps %>%
    pivot_wider(
      id_cols = c(gene, start, stop, strand, exonLength, CytCount),
      names_from = Context,
      values_from = ratio
    )

  # Convert to data.table and handle missing values
  wide <- as.data.table(wide)
  wide[, `NA` := NULL]
  wide$CG[is.na(wide$CG)] <- 0
  wide$CHH[is.na(wide$CHH)] <- 0
  wide$CHG[is.na(wide$CHG)] <- 0

  return(wide)
}
pctVarPlot <- function(PctVar) {

  # Clean the 'Model' column by removing '_ranked'
  PctVar$Model <- gsub("_ranked", "", PctVar$Model)

  # Arrange data by 'pct_var'
  PctVar <- PctVar %>% arrange(pct_var)

  # Define categories based on regex matching
  marks <- grep("H3|H2|H4|ATAC|CG|CHH|CHG|ATAC", PctVar$Model, value = TRUE)
  expn <- grep("avg|tissue|cv", PctVar$Model, value = TRUE)
  prot <- grep("IBM|SDG|Pol|JMJ", PctVar$Model, value = TRUE)

  # Set the factor levels for 'Model' based on categories
  PctVar$Model <- factor(PctVar$Model, levels = c(expn, prot, marks))

  # Create the plot
  PctVarMarks <- ggplot(PctVar, aes(x = -pct_var, y = Model)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.75) +
    theme_classic(base_size = 6) +
    labs(x = "% Var Explained") +
    scale_x_continuous(labels = abs, limits = c(-48, 4.5)) +
    scale_y_discrete(position = "right") +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(hjust = 0.5),
          axis.line.y = element_blank(),
          strip.text = element_blank())

  # Return the plot
  return(PctVarMarks)
}
corHeatmap <- function(DT) {

  # Create the correlation plot
  CorPlot <- ggplot(DT, aes(x = rn, y = mark, fill = value)) +
    geom_point(aes(size = `%IncMSE`), col = "black", shape = 21, stroke = 0.5) +
    coord_fixed() +
    theme_classic(base_size = 6) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(0.25, "cm"),
      #axis.text.y = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    scale_fill_gradient2(
      high = "#0570b0",
      low = "#cc4c02",
      mid = "#fff7fb",
      name = "Spearman \ncoefficient",
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
    )

  # Return the plot
  return(CorPlot)
}
intronGroupPlots <- function(Marks, epigmarks, column_name, ncol = 4) {
  # Dynamically set the column name to be used for "intronGroup"
  intron_column <- paste0(column_name)

  # Create the 'intronGroup' based on the specified column
  Marks$intronGroup <- cut2(Marks[[column_name]], g = 10)  # using the column_name dynamically

  # Get summaries (adjust this if needed for your function)
  all <- introngroup_summaries(Marks, epigmarks)

  # Generate the list of plots using lapply
  Plots <- lapply(epigmarks, function(mark) {
    plot_introngroup2(all, mark, intron_column)
  })

  # Combine plots into a grid
  col1 <- plot_grid(plotlist = Plots, ncol = ncol, align = "v")

  # Return the plot grid
  return(col1)
}

intronScatterPlots <- function(data, column, title, epigmarks) {
  # Ensure column is a string
  column <- as.character(column)

  # Extract the column as a numeric vector
  data$intronGroup <- cut2(as.numeric(data[[column]]), g = 10)

  # Step 2: Summarize data
  all <- introngroup_summaries(data, epigmarks) %>% ungroup() %>% as.data.frame()

  # Debugging step: Check the structure of `all`
  print(str(all))

  # Generate plots
  IntronNumber <- lapply(epigmarks, function(mark) plot_introngroup2(all, mark, title))

  return(IntronNumber)
}


#
#
#   # Step 3: Generate individual plots
#   IntronNumber <- lapply(epigmarks, function(mark) {
#     print(paste("Processing mark:", mark))  # Debugging check
#     plot_introngroup2(all, mark, title)
#   })
#
#   # Step 4: Adjust x-axis visibility for all but the last 4 plots
#   IntronNumber <- lapply(seq_along(IntronNumber), function(i) {
#     if (i <= length(IntronNumber) - 4) {
#       # Remove x-axis text and title for all but the last 4 plots
#       IntronNumber[[i]] + theme(axis.text.x = element_blank(),
#                                 axis.title.x = element_blank())
#     } else {
#       # Keep x-axis for the last 4 plots
#       IntronNumber[[i]] + theme(axis.title.x = element_blank())
#     }
#   })
#
#   # Step 5: Create title object
#   title_plot <- ggdraw() + draw_label(title, fontface = 'bold', x = 0, hjust = 0, size = 8)
#
#   # Step 6: Arrange plots
#   plot_grid_obj <- plot_grid(plotlist = IntronNumber,
#                              nrow = ceiling(length(IntronNumber) / 4),
#                              align = "v",
#                              rel_heights = c(rep(1, (length(IntronNumber) - 4) / 4), 1.25))
#
#   # Step 7: Combine title and plots
#   final_plot <- plot_grid(title_plot, plot_grid_obj, ncol = 1, align = "v", rel_heights = c(0.03, 1))
#
#   return(final_plot)
#}
plot_ttest_by_source <- function(src_value, Dups, MasterData, epigmarks) {
  all <- rbindlist(lapply(epigmarks, function(m) {
    result <- intron_marks_dups(Dups[src == src_value], MasterData, m)
    ttest <- result$ttest
    return(data.table(mark = m, statistic = ttest$statistic, p = ttest$p.value))
  }))

  all$mark <- factor(all$mark, levels = all$mark[order(all$statistic)])

  ttest_plot <- ggplot(all, aes(x = mark, y = -statistic, fill = -log10(p))) +
    geom_bar(stat = "identity", colour = "black") +
    theme_classic(base_size = 6) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      legend.key.size = unit(0.25, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    ) +
    scale_fill_gradient(
      low = "#edf8fb",
      high = "#0570b0",
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
    ) +
    labs(y = "T-test statistic")

  return(ttest_plot)
}
plot_ttest <- function( Dups, MasterData, epigmarks) {
  all <- rbindlist(lapply(epigmarks, function(m) {
    result <- intron_marks_dups(Dups, MasterData, m)
    ttest <- result$ttest
    return(data.table(mark = m, statistic = ttest$statistic, p = ttest$p.value))
  }))

  all$mark <- factor(all$mark, levels = all$mark[order(all$statistic)])

  ttest_plot <- ggplot(all, aes(x = mark, y = -statistic, fill = -log10(p))) +
    geom_bar(stat = "identity", colour = "black") +
    theme_classic(base_size = 6) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      legend.key.size = unit(0.25, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    ) +
    scale_fill_gradient(
      low = "#edf8fb",
      high = "#0570b0",
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
    ) +
    labs(y = "T-test statistic")

  return(ttest_plot)
}
plot_ttestTransposed <- function( Dups, MasterData, epigmarks) {
  all <- rbindlist(lapply(epigmarks, function(m) {
    result <- intron_marks_transposed(Dups, MasterData, m)
    ttest <- result$ttest
    return(data.table(mark = m, statistic = ttest$statistic, p = ttest$p.value))
  }))

  all$mark<-factor(all$mark, levels=all$mark[rev(order(all$statistic))])

  ttest_plot <- ggplot(all, aes(x=-statistic, y=mark, fill=-log10(p))) +
    geom_bar(stat="identity", colour="black")+
    theme_classic(base_size = 6)+
    theme(#axis.text.x=element_text(angle=45, hjust=1),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      legend.key.size = unit(0.25, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA))+
    scale_fill_gradient(low = "#edf8fb",high = "#cc4c02",
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    labs(x="T-test statistic")

  return(ttest_plot)
}
generate_correlation <- function(group_name) {
  # Create the group and data variables based on the passed group_name
  group <- column_groups[[group_name]]
  data <- ArabidopsisScores[, ..group]

  # Calculate Pearson correlation
  pearson_corr <- cor(data, method = "pearson", use = "pairwise.complete.obs")

  # Write the correlation to a CSV file with dynamic file name
  output_filename <- paste0(group_name, "_correlation.csv")
  fwrite(as.data.table(pearson_corr, keep.rownames = "Variable"),
         file = file.path("../Output-files/file-correlations", output_filename))
}
correlationDataTableWithPValues <- function(data_table, var1_columns, var2_columns) {
  # Initialize empty lists to store the correlation coefficients and p-values
  rho_list <- list()
  pvalue_list <- list()

  # Loop through each pair of var1_columns and var2_columns
  for (var1 in var1_columns) {
    for (var2 in var2_columns) {
      # Perform the correlation test (e.g., Pearson)
      #filteredData<-na.omit(data_table[, c(var1, var2)])
      cor_test_result <- cor.test(data_table[[var1]], data_table[[var2]], method = "spearman",exact = FALSE)

      # Append the correlation coefficient and p-value to the lists
      rho_list <- append(rho_list, list(cor_test_result$estimate))
      pvalue_list <- append(pvalue_list, list(cor_test_result$p.value))
    }
  }

  # Create two data tables, one for correlation coefficients and one for p-values
  rho_table <- data.table(Var1 = rep(var1_columns, each = length(var2_columns)),
                          Var2 = rep(var2_columns, times = length(var1_columns)),
                          Rho = unlist(rho_list))

  pvalue_table <- data.table(Var1 = rep(var1_columns, each = length(var2_columns)),
                             Var2 = rep(var2_columns, times = length(var1_columns)),
                             PValue = unlist(pvalue_list))

  return(list(RhoTable = rho_table, PValueTable = pvalue_table))
}
IntronGroupDF <- function(data, grouping_var, epigmarks) {
  # Group data
  data$intronGroup <- cut2(data[[grouping_var]], g = 10)

  # Summarize data
  all <- introngroup_summaries(data, epigmarks)
  all$Var<-grouping_var
  return(all)
}
IntronGroupSupplementaries <- function(data, grouping_var, label, epigmarks, ncol = 4) {
  # Group data
  data$intronGroup <- cut2(data[[grouping_var]], g = 10)

  # Summarize data
  all <- introngroup_summaries(data, epigmarks)

  # Generate all plots once
  plot_list <- lapply(epigmarks, function(mark) {
    ggplot(all[all$Mark == mark,], aes(x = intronGroup, y = Mean)) +
      geom_point(shape = 21, size = 1.8, fill = "#cc4c02") +
      geom_linerange(aes(ymin = Mean - SE, ymax = Mean + SE), color = "black", alpha = 0.6) +
      theme_classic(base_size = 6) +
      labs(x = label, y = mark) +
      theme(legend.position = "none")
  })

  # Version without x-axis text
  plot_list_no_x <- lapply(plot_list, function(p) {
    p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  })

  # Version with x-axis text only for last 4 plots
  plot_list_x <- lapply(seq_along(plot_list), function(i) {
    if (i <= length(plot_list) - 4) {
      plot_list[[i]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    } else {
      plot_list[[i]] + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })

  # Combine plots
  combined_plot_no_x <- plot_grid(plotlist = plot_list_no_x, ncol = ncol, align = "v")
  combined_plot_x <- plot_grid(plotlist = plot_list_x, ncol = ncol, align = "v")

  return(list(no_x = combined_plot_no_x, with_x = combined_plot_x, all=all))
}
