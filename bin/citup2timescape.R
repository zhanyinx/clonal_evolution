#!/usr/bin/env Rscript

library("getopt")
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "best_tree", "b", 1, "integer",
  "citup", "c", 1, "character",
  "mut_sample_clus", "m", 1, "character",
  "output", "o", 1, "character",
  "pyclone", "p", 1, "character",
  "sorting", "s", 1, "character",
  "help", "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}


library(rhdf5)
library(tidyverse)
library(reshape2)
library(timescape)
library(htmlwidgets)
library(stringr)

clone_freq <- h5read(opt$citup, paste0("/trees/",opt$best_tree,"/clone_freq"))
adjacency_list <- h5read(opt$citup,paste0("/trees/",opt$best_tree,"/adjacency_list"))
adj_list_df <- as.data.frame(t(as.data.frame.list(list(adjacency_list$block0_values))))
adj_list_df <- adj_list_df +1
colnames(adj_list_df) <- c("source","target")
clone_freq_df <- as.data.frame.list(list(clone_freq$block0_values))


tmp = read.delim(opt$mut_sample_clus, row.names=1)
colnames(clone_freq_df) <- colnames(tmp)

# sorting
if(file.exists(opt$sorting)){
  tmp = read.table(opt$sorting)
  tmp = tmp[order(tmp[,1]),]
  clone_freq_df <- clone_freq_df[,tmp$V2]
}else{
  warnings("samples are not sorted!")
}


mutations_tmp = read.delim(opt$pyclone)

conversion = h5read(opt$citup,paste0("/trees/",opt$best_tree,"/cluster_assignment"))
values = conversion$values + 1
names(values) = conversion$index

chrom = unlist(lapply(str_split(mutations_tmp$mutation_id, ":"), "[[", 1))
coord = unlist(lapply(str_split(unlist(
  lapply(str_split(mutations_tmp$mutation_id, ":"), "[[", 2)
), "-"), "[[", 1))
clone_id = values[as.character(mutations_tmp$cluster_id)]
timepoint = make.names(mutations_tmp$sample_id)
VAF = mutations_tmp$cellular_prevalence
genes = unlist(lapply(str_split(mutations_tmp$mutation_id, "\\."), "[[", 4))
genes = paste0(genes, VAF, sep="_")
mutations = data.frame(chrom=chrom, coord=coord, clone_id=clone_id, timepoint=timepoint, VAF=VAF, genes = genes)


clone_prev_general <- setNames(melt((as.matrix(clone_freq_df)), value.name="clonal_prev"), c("clone_id","timepoint", "clonal_prev"))
reshaped_clone_prevs <- clone_prev_general[,c(2,1,3)]
obj = timescape(clonal_prev=reshaped_clone_prevs, tree_edges=adj_list_df, mutations=mutations)
# obj = timescape(clonal_prev=reshaped_clone_prevs, tree_edges=adj_list_df)

saveWidget(obj, file=paste0(opt$output,".html"))
