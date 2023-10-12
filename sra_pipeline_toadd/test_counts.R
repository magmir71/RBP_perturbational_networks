# edgeR count normalization and natural breaks clustering #
library(Ckmeans.1d.dp, quietly = T)
library(RColorBrewer, quietly = T)
library(vroom, quietly = T)
library(tidyverse, quietly = T)
library(edgeR, quietly = T)
library(DESeq2, quietly = T)

#### INPUT ####
args <- commandArgs(trailingOnly = T)

outdir <- args[1]
sample_data <- args[2]
genes.bed <- args[3]
sample_dir <- args[4]
path_to_GTF_auc <- 'BAMS/withGTF/megadepth/genes'
path_to_noGTF_auc <- 'BAMS/noGTF/megadepth/genes'

#########
#########

sample_data <- read.delim(sample_data, header = T, sep = ',')
sample_data$LibraryLayout[sample_data$LibraryLayout == 'PAIRED'] <- 'Paired'
sample_data$LibraryLayout[sample_data$LibraryLayout == 'SINGLE'] <- 'Single'
genes.bed <- read.delim(genes.bed, header = F)

# add gene count files paths_withGTF as column
sample_data$paths_withGTF <- file.path(sample_dir, sample_data$LibraryLayout, sample_data$sample, path_to_GTF_auc, paste0(sample_data$sample, '_auc.out'))
# prepare list for reading
list_df_withGTF <- with(sample_data, split(paths_withGTF, experiment_id))


#### GET GENE IDS FROM GTF FOR COUNTS FILES####
x <- head(read.delim(sample_data$paths_withGTF[1], header = F, sep = '\t'), -2)[,-4]
genenames <- merge(x, genes.bed[,1:4], by = c('V1', 'V2', 'V3'), sort = F)
genenames <- genenames[!duplicated(genenames),] #for some reason there's a couple of dup regions
colnames(genenames) <- c('chr', 'start', 'end', 'gene_id')
genenames$gene_id <- gsub("\\.[0-9]+","", genenames$gene_id)

#############

#### MEGADEPTH COUNTS FILES READING ####

# creation of temp list of runs for each experiment
for (i in sample_data$experiment_id) {
  j <- as.list(unlist(as.list(list_df_withGTF[i])))
  names(j) <- sample_data$sample[sample_data$experiment_id==i]
  assign(i, list(j))
}

# aggregation of those lists into one single list (one for GTF and another for noGTF (later))
final_list_withGTF <- list()
for (i in ls(pattern='SRP')) {
  final_list_withGTF <- append(final_list_withGTF, get(i))
}
names(final_list_withGTF) <- unique(sample_data$experiment_id)

# rm of temp lists
for (i in ls(pattern='SRP')) {
  rm(i)
}

#### MERGING PER EXPERIMENT - CPM NORM - TARGET GENE PLOT ####
# exp design
design <- sample_data[,c('sample','EXP_CTL')]

setwd(outdir)

# loop in final_list_withGTF for merge - norm - plot
# with gtf
for (i in final_list_withGTF) {
  merged_df <- as.data.frame(map_dfc(i,
                                     ~vroom(.x,
                                            col_select=c(4),
                                            delim = '\t',
                                            col_names = FALSE,
                                            show_col_types = FALSE)))
  #renaming cols
  suppressWarnings(colnames(merged_df) <- names(i))
  suppressWarnings(merged_df <- head(merged_df, -2)) #remove last two lines
  rownames(merged_df) <- genenames$gene_id #gene ids as rownames
  
  #cat(unique(sample_data$experiment_id[sample_data$sample==names(i)]))
  cat('Merging of raw counts within experiments...\n')
  #cat(head(merged_df))
  
  cat('CPM normalization within experiments...\n')
  cat(calcNormFactors(merged_df),'\n')
  cpm_counts <- cpm(merged_df)
  #cat(head(cpm_counts))
  
  target_counts <- cpm_counts[rownames(cpm_counts)==unique(sample_data$targeted_gene[sample_data$sample==names(i)])]
  cat('Target gene CPM:\n')
  cat(unique(sample_data$targeted_gene[sample_data$sample==names(i)]), '\n')
  cat(target_counts)
  
  k <- 2
  result <- Ckmeans.1d.dp(target_counts, k)
  colors <- brewer.pal(k, "Dark2")
  midpoints <- ahist(result, style="midpoints", data=target_counts, plot=FALSE)$breaks[2:k]
  
  cat('\nCreatin pdf...\n')
  pdf(paste0('Clustering_plot_', unique(sample_data$experiment_id[sample_data$sample==names(i)]), '_GTF.pdf'), height = 2.5)
  p <- plot(result, 
            col.clusters = colors, 
            main = paste('K-means clustering of', unique(unique(sample_data$experiment_id[sample_data$paths_withGTF==i]))),
            ylab = '',
            yaxt = 'n',
            xlab = paste(unique(sample_data$targeted_gene[sample_data$sample==names(i)]), 'CPM', '(GTF)'),
            ylim=c(0,1))
  abline(v=midpoints, col="RoyalBlue", lwd=3, lty=3,)
  text(x = target_counts, y = par("usr")[3]-0.15, labels = paste(names(i), design$EXP_CTL[design$sample==names(i)]), xpd = NA, srt=45, cex=0.4, adj = 0.9)
  dev.off()
}


#### NOGTF ####

# add gene count files paths_withGTF as column
sample_data$paths_noGTF <- file.path(sample_dir, sample_data$LibraryLayout, sample_data$sample, path_to_noGTF_auc, paste0(sample_data$sample, '_auc.out'))
# prepare list for reading
list_df_noGTF <- with(sample_data, split(paths_noGTF, experiment_id))

for (i in sample_data$experiment_id) {
  j <- as.list(unlist(as.list(list_df_noGTF[i])))
  names(j) <- sample_data$sample[sample_data$experiment_id==i]
  assign(i, list(j))
}

final_list_noGTF <- list()
for (i in ls(pattern='SRP')) {
  final_list_noGTF <- append(final_list_noGTF, get(i))
}
names(final_list_noGTF) <- unique(sample_data$experiment_id)

#### no GTF loop ####
for (i in final_list_noGTF) {
  merged_df <- as.data.frame(map_dfc(i,
                                     ~vroom(.x,
                                            col_select=c(4),
                                            delim = '\t',
                                            col_names = FALSE,
                                            show_col_types = FALSE)))
  #renaming cols
  suppressWarnings(colnames(merged_df) <- names(i))
  suppressWarnings(merged_df <- head(merged_df, -2)) #remove last two lines
  rownames(merged_df) <- genenames$gene_id #gene ids as rownames
  
  #cat(unique(sample_data$experiment_id[sample_data$sample==names(i)]))
  cat('Merging of raw counts within experiments...\n')
  #cat(head(merged_df))
  
  cat('CPM normalization within experiments...\n')
    cat(calcNormFactors(merged_df), '\n')
  cpm_counts <- cpm(merged_df)
  cat('\n')
  #cat(head(cpm_counts))
  
  target_counts <- cpm_counts[rownames(cpm_counts)==unique(sample_data$targeted_gene[sample_data$sample==names(i)])]
  cat('Target gene CPM:\n')
  cat(unique(sample_data$targeted_gene[sample_data$sample==names(i)]), '\n')
  cat(target_counts)

  k <- 2
  result <- Ckmeans.1d.dp(target_counts, k)
  colors <- brewer.pal(k, "Dark2")
  midpoints <- ahist(result, style="midpoints", data=target_counts, plot=FALSE)$breaks[2:k]
  
  cat('\nCreatin pdf...\n')
  pdf(paste0('Clustering_plot_', unique(sample_data$experiment_id[sample_data$sample==names(i)]), '_noGTF.pdf'), height = 2.5)
  p <- plot(result, 
            col.clusters = colors, 
            main = paste('K-means clustering of', unique(unique(sample_data$experiment_id[sample_data$paths_noGTF==i]))),
            ylab = '',
            yaxt = 'n',
            xlab = paste(unique(sample_data$targeted_gene[sample_data$sample==names(i)]), 'CPM', '(noGTF)'),
            ylim=c(0,1))
  abline(v=midpoints, col="RoyalBlue", lwd=3, lty=3,)
  text(x = target_counts, y = par("usr")[3]-0.15, labels = paste(names(i), design$EXP_CTL[design$sample==names(i)]), xpd = NA, srt=45, cex=0.4, adj = 0.9)
  dev.off()
}



#first cluster
#width
#max(target_counts[result$cluster==1]) - min(target_counts[result$cluster==1])
#
#distance from midpoint
#midpoints - max(target_counts[result$cluster==1])
#
#second cluster
#width
#max(target_counts[result$cluster==2]) - min(target_counts[result$cluster==2])
#
#distance from midpoint
#midpoints - min(target_counts[result$cluster==2])

