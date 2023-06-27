#!/usr/bin/env Rscript

library('getopt')
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose'            , 'v', 2, "integer",
  'help'               , 'h', 0, "logical",
  'pyclone_loci'       , 'p', 1, "character",
  'output_directory'   , 'o', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# Check if pyclone_loci file exists
if (!file.exists(opt$pyclone_loci)) {
  stop("Error: pyclone_loci file does not exist.")
}

# Check if output_directory exists and is a directory
if(!dir.exists(opt$output_directory)){
  # Create the output_directory if it doesn't exist
  dir.create(opt$output_directory, recursive = TRUE)
}


library(tidyverse)
library(reshape2)

pyclone2citup <- function(pyclone_loci,output) { 

        l <- read_tsv(pyclone_loci)

        write_tsv(as.data.frame(max(l$cluster_id)), paste(output,"num_clusters.tsv",sep = "/"),col_names =F)

        casted.l <- acast(l, mutation_id~sample_id, value.var='cellular_prevalence')
        casted.clust_mut <- acast(l, cluster_id~sample_id, value.var='cellular_prevalence')
        #View(casted.clust_mut)
        write.table(casted.clust_mut, file = paste(output,"num_mut_clust_id.txt", sep = "/"),sep = "\t",quote = F)

        l_tib <- as_tibble(casted.l)
        write_tsv(l_tib,paste(output,"citup_frequencies.tsv", sep = "/"), col_names= F)
        casted.clus<-(acast(l, mutation_id~sample_id, value.var='cluster_id'))
        write.table(casted.clus, file = paste(output,"mut_sample_clus.txt", sep = "/"),sep = "\t",quote = F)

        casted.mut_clus<-(acast(l, mutation_id~cluster_id, value.var='cellular_prevalence'))
        write.table(casted.mut_clus, file = paste(output,"mut_id_clust_id.txt",sep = "/"),quote = F, sep = "\t")

        write_tsv((as.data.frame(casted.clus))[1],paste(output,"citup_clusters.tsv",sep = "/"),col_names = F)
        write_tsv((as.data.frame(casted.clus))[1]-as.integer(1),paste(output,"citup_clusters_0based.tsv",sep = "/"),col_names = F)

}

pyclone2citup(opt$pyclone_loci, opt$output_directory)