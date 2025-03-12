#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
suppressPackageStartupMessages(library(zoo))
library(ggpubr)
suppressPackageStartupMessages(library(tidyverse))
library(dplyr, exclude = c("filter", "lag"))
library(tidyr)
library(optparse)


# Define options
option_list <- list(
make_option(c("-x", "--refstats1"), type = "character",
help = "Refstats (mean ionic current) file for sample 1 [required]"),
make_option(c("-y", "--refstats2"), type = "character",
help = "Refstats (mean ionic current) file for sample 2 [required]"),
make_option(c("-k", "--ksinput"), type = "character",
help = "Refstats KS test file [required]"),
make_option(c("-p", "--pvalue"), type = "numeric", default = "0.05",
help = "Cutoff for generalized ESD p-value [default 0.05]"),
make_option(c("-c", "--ksnorm"), type = "numeric", default = "3",
help = "Cutoff value for the normalized KS statistic [default 3]"),
make_option(c("-n", "--name"), type = "character", default = "DRAMA",
help = "Chromosome name or file prefix [default %default]"),
make_option(c("-l","--labels"), action="store_true", default=FALSE,
help="Include top 10 outlier labels on the plot (default is FALSE)"),
make_option(c("-d","--height"), type="numeric", default="5",
help="User-provided height for each plot (default=5)"),
make_option(c("-w","--width"), type="numeric", default="6",
help="User-provided width for each plot (default=6)")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
options <- parse_args(opt_parser)
add_labels <- options$labels

# Load files
cat("\n\nLoading input files...\n\n\n")
#sample_1 <- read.delim(args[1], sep="\t", header=TRUE)
#sample_2 <- read.delim(args[2], sep="\t", header=TRUE)
#KS <- read.delim(args[3], sep="\t", header=TRUE)
#chrom=as.character(args[4])

sample_1 <- read.delim(options$refstats1, sep="\t", header=TRUE)
sample_2 <- read.delim(options$refstats2, sep="\t", header=TRUE)
KS <- read.delim(options$ksinput, sep="\t", header=TRUE)
chrom=as.character(options$name)

# separate positive/negative strand
sample_1_pos <- sample_1[(sample_1$strand == "+"),]
sample_2_pos <- sample_2[(sample_2$strand == "+"),]
sample_1_neg <- sample_1[(sample_1$strand == "-"),]
sample_2_neg <- sample_2[(sample_2$strand == "-"),]
KS_pos <- KS[(KS$strand == "+"),]
KS_neg <- KS[(KS$strand == "-"),]

# Check if both strands are present in KS file
if(length(KS_pos$ref) < 1) {
    cat("\n*** No positive strand values in refstats file, output will be negative strand only ***\n\n")
}

if(length(KS_neg$ref) < 1) {
    cat("\n*** No negative strand values in refstats file, output will be positive strand only ***\n\n")
}


# Set function: generate an index column for breaking into transcript regions

indexRegions <- function(combined_df) {
    colnames(combined_df) <- c("chrom","ref","sample1_depth", # rename columns
                                "sample1_current_mean",
                                "sample2_depth",
                                "sample2_current_mean",
                                "ks_stat")
    combined_df$group <- c(1, cumsum(diff(combined_df$ref) >= 25) + 1)  # Create a group ID based on consecutive indices
    return(combined_df)
}


# Set function: loop through each group (transcript region) in index and normalize KS stat and mean ionic current for depth

normalizeStats <- function(indexed_df) {
    window_size <- 65
    central_position <- floor(window_size/2) +1
    ## KS statistic ##
    exclude_n <- 10 # Set number of bases on either end of each region/gene to exclude due to noise from a partial sliding window
    indexed_df$ks_stat.window.means <- NA # Create a new column to store rolling mean for KS stat
    for(i in unique(indexed_df$group)) {
        group_data <- indexed_df[indexed_df$group == i, ] # Subset the data for the current group
        if (nrow(group_data) > 2 * exclude_n) { # make sure group size is large enough to remove 10 bases
            group_data_to_apply <- group_data[(exclude_n + 1):(nrow(group_data) - exclude_n), ] # exclude bases on either end
            values <- rollapply( # calculate sliding window means
                group_data_to_apply$ks_stat.prior.count,
                width = window_size,
                FUN = function(x) {
                    x_remove_center <- x[-central_position]
                    mean(x_remove_center)
                },
                align = "center",
                fill = NA, # fill with NA if position is excluded from analysis (10 bases on either end of window)
                partial = TRUE
            )
            # Update the main dataframe with the results
            indexed_df$ks_stat.window.means[indexed_df$group == i][(exclude_n + 1):(nrow(group_data) - exclude_n)]  <- values
        }
    }
    indexed_df <- indexed_df[!is.na(indexed_df$ks_stat.window.means), ] # remove any "NA" values
    indexed_df$ks_stat.mean.window.FC <- indexed_df$ks_stat.prior.count / indexed_df$ks_stat.window.means # divide the KS statistic by the 65-base window mean
    indexed_df <- indexed_df[!is.na(indexed_df$ks_stat.mean.window.FC),] # remove any "NA" values

    ## Mean ionic current ##
    indexed_df$current.mean.diff <- indexed_df$sample1_current_mean - indexed_df$sample2_current_mean # get the difference between the mean ionic current for both samples
    
    cat(paste0("Calculating Generalized ESD for mean ionic current with alpha=",options$pvalue,"...\n\n"))
    
    for(i in unique(indexed_df$group)) {
        group_data <- indexed_df[indexed_df$group == i, ] # Subset the data for the current group
        if (nrow(group_data) > 2 * exclude_n) { # make sure group size is large enough to remove 10 bases
            values <- rollapply( # calculate sliding window means
                group_data$current.mean.diff,
                width = window_size,
                FUN = function(x) {
                    x_remove_center <- x[-central_position]
                    mean(abs(x_remove_center))
                },
                align = "center",
                fill = NA, # fill with NA if position is excluded from analysis (10 bases on either end of window)
                partial = TRUE
            )
            values[1:exclude_n] <- NA
            values[(length(values) - exclude_n + 1):length(values)] <- NA
            # Update the main dataframe with the results
            indexed_df$current.window.means.FC[indexed_df$group == i] <-
                log2(abs(group_data$current.mean.diff) / values)

            # Calculate generalized ESD for values in index group
            subset_range <- (exclude_n + 1):(nrow(group_data) - exclude_n)
            norm_current <- indexed_df$current.window.means.FC[indexed_df$group == i][subset_range]
            alpha <- options$pvalue
            outliers <- numeric(0)
            len <- length(norm_current)
            mean_data <- mean(norm_current)
            sd_data <- sd(norm_current)
            p_values <- numeric(len) # Vector to store p-values
            right_tailed_flags <- logical(len) # Vector to store right-tailed flags (outliers)
            right_tailed_all <- logical(len) # Vector to store all right-tailed points (anything above the mean)
            for (n in 1:len) {
                residual <- abs(norm_current[n] - mean_data) / sd_data # Calculate the studentized residual for each data point
                # Calculate the p-value based on the residual
                # Use the t-distribution with (len - 1) degrees of freedom
                p_value <- 2 * (1 - pt(residual, df = len - 1))  # Two-tailed p-value
                # Store the p-value
                p_values[n] <- p_value
                # Flag as a right-tailed outlier if p-value is less than alpha (one-tailed test)
                if (!is.na(p_value) && p_value < alpha && norm_current[n] > mean_data) {    
                    right_tailed_flags[n] <- TRUE
                } else {
                    right_tailed_flags[n] <- FALSE
                }
                # Flag as right-tailed (not just outliers) for plotting purposes
                if (norm_current[n] > mean_data) {
                    right_tailed_all[n] <- TRUE
                } else {
                    right_tailed_all[n] <- FALSE
                }
            }
            # Update p-values and flags
            indexed_df$esd.pval[indexed_df$group == i] <- NA
            indexed_df$right_tailed[indexed_df$group == i] <- NA
            indexed_df$right_tailed_outlier[indexed_df$group == i] <- NA
            indexed_df$esd.pval[indexed_df$group == i][(exclude_n + 1):(nrow(group_data) - exclude_n)] <- p_values
            indexed_df$right_tailed[indexed_df$group == i][(exclude_n + 1):(nrow(group_data) - exclude_n)] <- right_tailed_all
            indexed_df$right_tailed_outlier[indexed_df$group == i][(exclude_n + 1):(nrow(group_data) - exclude_n)] <- right_tailed_flags
        } else {
            indexed_df$current.window.means.FC[indexed_df$group == i] <- NA
            indexed_df$esd.pval[indexed_df$group == i] <- NA
            indexed_df$right_tailed[indexed_df$group == i] <- NA
            indexed_df$right_tailed_outlier[indexed_df$group == i] <- NA
        }
    }
    indexed_df <- na.omit(indexed_df)
    
    return(indexed_df)
}

# Set function: collapse rows with consecutive reference positions

collapseConsecutive <- function(diff_mod_df) {
    diff_mod_df <- diff_mod_df %>% 
    arrange(ref) %>%
    mutate(
        #ref_block = cumsum(ref != "ref" | dplyr::lag(ref != "ref", default = TRUE))
        ref_block = cumsum(ref != lag(ref, default = first(ref)) + 1)
        )
    diff_mod_df_collapsed <- diff_mod_df %>%
        group_by(ref_block) %>%
        mutate(
            esd.pval_min = ifelse(ref == "ref", min(esd.pval), esd.pval)
            ) %>%
            dplyr::filter(ref != "ref" | row_number() == 1) %>%
            distinct(ref_block, .keep_all = TRUE) %>%
            ungroup()
    return(diff_mod_df_collapsed)
}

options(scipen=999) # make sure format of bedgraph files will work with IGV

# Positive strand analysis
if(length(KS_pos$ref) > 0) {
    cat("\n--------------------------------------\n
Performing positive strand analysis...\n
--------------------------------------\n\n")
    merge_stats_pos <- merge(sample_1_pos,sample_2_pos, by='ref') # merge mean current stats from both samples
    combined_all_col <- merge(merge_stats_pos,KS_pos, by='ref') # merge with KS stat into a single dataframe
    combined_pos <- combined_all_col[,c(2,1,4,5,8,9,14)] # filter out unnecessary columns
    cat("Indexing transcript regions...\n\n")
    combined_pos <- indexRegions(combined_pos)
    combined_pos$ks_stat.prior.count <- combined_pos$ks_stat + 0.01 # add prior count to KS statistic to avoid effects of zero values
    cat("Normalizing mean ionic current and KS statistic using sliding windows...\n\n")
    combined_pos <- normalizeStats(combined_pos)
    combined_pos_righttail <- combined_pos[(combined_pos$right_tailed==TRUE),]
    combined_pos_righttail <- as.data.frame(combined_pos_righttail)
    combined_pos$ref2 <- combined_pos$ref + 1
    combined_pos_righttail$ref2 <- combined_pos_righttail$ref + 1
    # Filter with user-provided p-value and normalized KS statistic cutoffs
    cat(paste0("Calling differential modifications based on user-provided cutoff values:\n
           p-value < ",options$pvalue,"\n
           Normalized KS statistic > ",options$ksnorm,"\n\n"))
    diffmods_pos <- combined_pos[((combined_pos$right_tailed_outlier==TRUE)&(combined_pos$ks_stat.mean.window.FC>options$ksnorm)),]
    diffmods_pos <- collapseConsecutive(diffmods_pos)
    #right_tails_pos <- combined_pos[(combined_pos$right_tailed==TRUE),]
    combined_pos_righttail$diff_mod <- "no" # create a column for differential modification call
    combined_pos_righttail$diff_mod[((combined_pos_righttail$right_tailed_outlier==TRUE)&(combined_pos_righttail$ks_stat.mean.window.FC>3))] <- "yes" # change differential modifications to "yes"
    
    # Generate volcano plot and bed files for positive strand
    plot_filename=paste0(chrom, "_diffmods_pos.jpeg")
    plot_filename2=paste0(chrom, "_diffmods_pos_density.jpeg")
    bed_filename=paste0(chrom, "_diffmods_pos.bed", sep="")
    if (nrow(diffmods_pos) > 0) {
        cat("Plotting positive strand results...\n\n")
        top_10_extreme <- combined_pos_righttail %>%
            dplyr::filter(esd.pval < options$pvalue & ks_stat.mean.window.FC > options$ksnorm) %>%
            arrange(esd.pval) %>%
            slice_head(n = 25) %>%
            arrange(desc(ks_stat.mean.window.FC)) %>%
            slice_head(n = 10)  # Get the top 10 extreme points
        plot <- ggplot(combined_pos_righttail, aes(y=-log10(esd.pval),
                                        x=ks_stat.mean.window.FC,                
                                        col=diff_mod)) +
             geom_vline(xintercept = options$ksnorm, col = "gray", linetype = 'dashed') +
             geom_hline(yintercept = -log10(options$pvalue), col = "gray", linetype = 'dashed') +
             geom_point(size = 1) + 
             theme_minimal() +
             scale_color_manual(values = c("#B2BEB5", "#bb0c00"),
                                labels = c("Not significant", "Differentially\nmodified")) +
             xlab("\nNormalized KS statistic") +
             ylab("-Log10 p-value: normalized\nmean ionic current difference\n") +
             theme(legend.position="none",
                    plot.background = element_rect(
                    fill = "white",
                    colour = "white")) +
             xlim(0,max(combined_pos_righttail$ks_stat.mean.window.FC)) +
             ylim(0,max(-log10(combined_pos_righttail$esd.pval)))
        if (add_labels == TRUE) {
            plot <- plot + 
            geom_text_repel(data=top_10_extreme,
                        aes(label=ref2, fontface="bold"),size=3, show.legend = FALSE)
        }
        ggsave(plot_filename,width = options$width, plot=plot, height = options$height, dpi = 300)
                ggplot(combined_pos_righttail, aes(y=-log10(esd.pval),
                                        x=ks_stat.mean.window.FC)) +
                geom_vline(xintercept = options$ksnorm, col = "gray", linetype = 'dashed') +
                geom_hline(yintercept = -log10(options$pvalue), col = "gray", linetype = 'dashed') +
                theme_minimal() +
                xlab("\nNormalized KS statistic") +
                ylab("-Log10 p-value: normalized\nmean ]ionic current difference\n") +
                theme(plot.background = element_rect(
                    fill = "white",
                    colour = "white")) +
                scale_fill_gradient(low = "#2B2D85",high = "#FFEF1D") +
                xlim(0,max(combined_pos_righttail$ks_stat.mean.window.FC)) +
                ylim(0,max(-log10(combined_pos_righttail$esd.pval))) +
                geom_hex(bins=75, na.rm=TRUE)
        ggsave(plot_filename2,width = options$width, height = options$height, dpi = 300)
        mods_bed <- diffmods_pos[,c(1,2,17)] # format peaks table for bed file output
        mods_bed$name <- "diff_modification"
        mods_bed$pval <- diffmods_pos$esd.pval
        mods_bed$strand <- "+"
        mods_bed$ks_norm <- diffmods_pos$ks_stat.mean.window.FC
        write.table(mods_bed, bed_filename, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
        } else {
            cat("\n*** No differential modifications detected for positive strand.***\n\n\n")
        }
    all_bed_filename=paste0(chrom, "_all_pos.bed", sep="")
    all_bed <- combined_pos[,c(1,2,17)]
    all_bed$name <- "DRAMA_stats"
    all_bed$pval <- combined_pos$esd.pval
    all_bed$strand <- "+"
    all_bed$ks_norm <- combined_pos$ks_stat.mean.window.FC
    write.table(all_bed, all_bed_filename, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
}



# Negative strand analysis

if(length(KS_neg$ref) > 0) {
    cat("\n--------------------------------------\n
Performing negative strand analysis...\n
--------------------------------------\n\n")
    merge_stats_neg <- merge(sample_1_neg,sample_2_neg, by='ref') # merge mean current stats from both samples
    combined_all_col <- merge(merge_stats_neg,KS_neg, by='ref') # merge with KS stat into a single dataframe]
    combined_neg <- combined_all_col[,c(2,1,4,5,8,9,14)] # filter out unnecessary columns
    cat("Indexing transcript regions...\n\n")
    combined_neg <- indexRegions(combined_neg)
    combined_neg$ks_stat.prior.count <- combined_neg$ks_stat + 0.01 # add prior count to KS statistic to avoid effects of zero values
    cat("Normalizing mean ionic current and KS statistic using sliding windows...\n\n")
    combined_neg <- normalizeStats(combined_neg)
    combined_neg_righttail <- combined_neg[(combined_neg$right_tailed==TRUE),]
    combined_neg_righttail <- as.data.frame(combined_neg_righttail)
    combined_neg$ref2 <- combined_neg$ref + 1
    combined_neg_righttail$ref2 <- combined_neg_righttail$ref + 1
    # Filter with user-provided p-value and normalized KS statistic cutoffs
    cat(paste0("Calling differential modifications based on user-provided cutoff values:\n
           p-value < ",options$pvalue,"\n
           Normalized KS statistic > ",options$ksnorm,"\n\n"))
    diffmods_neg <- combined_neg[((combined_neg$right_tailed_outlier==TRUE)&(combined_neg$ks_stat.mean.window.FC>options$ksnorm)),]
    diffmods_neg <- collapseConsecutive(diffmods_neg)
    combined_neg_righttail$diff_mod <- "no" # create a column for differential modification call
    combined_neg_righttail$diff_mod[((combined_neg_righttail$right_tailed_outlier==TRUE)&(combined_neg_righttail$ks_stat.mean.window.FC>3))] <- "yes" # change differential modifications to "yes"

    # Generate volcano plot and bed files for negative strand
    plot_filename=paste0(chrom, "_diffmods_neg.jpeg")
    plot_filename2=paste0(chrom, "_diffmods_neg_density.jpeg")
    bed_filename=paste0(chrom, "_diffmods_neg.bed", sep="")
    if (nrow(diffmods_neg) > 0) {
        cat("Plotting negative strand results...\n\n")
        top_10_extreme_neg <- combined_neg_righttail %>%
            dplyr::filter(esd.pval < options$pvalue & ks_stat.mean.window.FC > options$ksnorm) %>%
            arrange(esd.pval) %>%
            slice_head(n = 25) %>%
            arrange(desc(ks_stat.mean.window.FC)) %>%
            slice_head(n = 10)  # Get the top 10 extreme points
        plot <- ggplot(combined_neg_righttail, aes(y=-log10(esd.pval),
                                        x=ks_stat.mean.window.FC,
                                        col=diff_mod)) +
             geom_vline(xintercept = options$ksnorm, col = "gray", linetype = 'dashed') +
             geom_hline(yintercept = -log10(options$pvalue), col = "gray", linetype = 'dashed') +
             geom_point(size = 1) +
             theme_minimal() +
             scale_color_manual(values = c("#B2BEB5", "#bb0c00")) +
                                #labels = c("Not significant", "Differentially\nmodified")) +
             xlab("\nNormalized KS statistic") +
             ylab("-Log10 p-value: normalized\nmean ionic current difference\n") +
             theme(legend.position="none",
                    plot.background = element_rect(
                    fill = "white",
                    colour = "white"))
        if (add_labels == TRUE) {
            plot <- plot +
            geom_text_repel(data=top_10_extreme,
                        aes(label=ref2, fontface="bold"),size=3, show.legend = FALSE)
        }
        ggsave(plot_filename,width = 6, plot=plot, height = 5, dpi = 300)
        ggplot(combined_neg_righttail, aes(y=-log10(esd.pval),
                                        x=ks_stat.mean.window.FC)) +
                geom_vline(xintercept = options$ksnorm, col = "gray", linetype = 'dashed') +
                geom_hline(yintercept = -log10(options$pvalue), col = "gray", linetype = 'dashed') +
                theme_minimal() +
                xlab("\nNormalized KS statistic") +
                ylab("-Log10 p-value: normalized\nmean ]ionic current difference\n") +
                theme(plot.background = element_rect(
                    fill = "white",
                    colour = "white")) +
                scale_fill_gradient(low = "#2B2D85",high = "#FFEF1D") +
                xlim(0,max(combined_neg_righttail$ks_stat.mean.window.FC)) +
                ylim(0,max(-log10(combined_neg_righttail$esd.pval))) +
                geom_hex(bins=75, na.rm=TRUE)
        ggsave(plot_filename2,width = 6.5, height = 5, dpi = 300)
        mods_bed <- diffmods_neg[,c(1,2,17)] # format peaks table for bed file output
        mods_bed$name <- "diff_modification"
        mods_bed$pval <- diffmods_neg$esd.pval
        mods_bed$strand <- "+"
        mods_bed$ks_norm <- diffmods_neg$ks_stat.mean.window.FC
        write.table(mods_bed, bed_filename, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
        } else {
            cat("\n***No differential modifications detected for negative strand.***\n\n\n")
        }
    all_bed_filename=paste0(chrom, "_all_neg.bed", sep="")
    all_bed <- combined_neg[,c(1,2,17)]
    all_bed$name <- "DRAMA_stats"
    all_bed$pval <- combined_neg$esd.pval
    all_bed$strand <- "+"
    all_bed$ks_norm <- combined_neg$ks_stat.mean.window.FC
    write.table(all_bed, all_bed_filename, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
}
