#!/usr/bin/env R

# USAGE: assemble_raw_results.R 
# Looks for .RData files from parallel simulations in current working directory

# Join results from parallel simulation of null distribution of
# proband ascertainment bias in GO enrichment.

# Builds table that is used to plot mean and standard deviation
# of q-values expected by random chance in stratification by
# PRS approach

rm(list = ls())

# load libraries
library(plyr)

# glob in RData files
filenames <- Sys.glob("*.RData")

# empty sapce to store master resutls list and iteration ID
master_results = list()
iterations = c()
master_samples = list()
master_genes = list()
k = -1 # for indexing into lists
f = 0 # for indexing filenames

# unpack RData objects, 100 total each with a unique filename
for(i in filenames){ 
    # increment counters
    f = f + 1
    k = k + 1

    # load file into  memory
    load(i)

    # loop through the 100 simulations in each of the files
    for(j in 1:100){
        # add file name to list of iterations
        iterations = c(iterations, paste(strsplit(filenames[f], split = "_")[[1]][1], j, sep = "_"))
        
        # append to master list
        master_results[[100*k+j]] = mc_results[[j]]

        # add samples in the iteration to the master sample list
        master_samples[[100*k+j]] = as.character(samples[[j]])
    
        # get genes for each iteration and save to master
        master_genes[[100*k+j]] = as.character(genes[[j]])
    }
}

# check length 10000
if(length(master_results) != 10000){
    stop("incorrect number of simulations loaded")
}

# get list of GO terms
uniq_go = master_results[[1]][,1]

# print number of go terms included in analysis
length(uniq_go)

# SAMPLE TABLE
# make table of iteration and sample IDs
sample_dict = as.data.frame(matrix(nrow = 325, ncol = length(iterations)))
colnames(sample_dict) = iterations
for(i in 1:length(iterations)){
    sample_dict[,i] = master_samples[[i]]
}

# save to file
write.table(sample_dict, file = "topGO_10K_sample_dict.txt", quote = FALSE, sep = "\t")

# GENE TABLE
# make table of genes in each simulation
for(i in 1:length(master_genes)){
    if(i == 1){
        gene_table = t(as.data.frame(master_genes[[i]]))
        colnames(gene_table) = gene_table[1,]
    } else {
        x = t(as.data.frame(master_genes[[i]]))
        colnames(x) = x[1,]
        gene_table = rbind.fill(as.data.frame(gene_table), as.data.frame(x))
    }
}

# set row names (stripped by plyr)
row.names(gene_table) = iterations

# convert to numeric (columns are factors)
for(i in 1:ncol(gene_table)){
    gene_table[,i] = as.numeric(gene_table[,i])
}

# convert NA to 0
gene_table[is.na(gene_table)] = 0

# write to file
write.table(gene_table, file = "topGO_10K_gene_table.txt", quote = FALSE, sep = "\t")

# P AND Q VALUES
# empty sapce for raw p-values: rows are iterations and columns are GO terms
df_pval = as.data.frame(matrix(nrow = length(master_results), ncol = length(uniq_go)))
colnames(df_pval) = uniq_go
row.names(df_pval) = iterations

# empty sapce for raw q-values: rows are iterations and columns are GO terms
df_qval = as.data.frame(matrix(nrow = length(master_results), ncol = length(uniq_go)))
colnames(df_qval) = uniq_go
row.names(df_qval) = iterations

# loop through each element in the master list, extracting the p and q values
# in the correct order and assigning as rows in new data frames
for(i in 1:length(master_results)){
    tmp = master_results[[i]]
    tmp = tmp[match(uniq_go, tmp$GO.ID) ,]
          
    df_pval[i,] = as.numeric(tmp$classic)
    df_qval[i,] = as.numeric(tmp$qval)
}

# save to file
write.table(df_pval, file = "topGO_10K_raw_pvals.txt", quote = FALSE, sep = "\t")
write.table(df_qval, file = "topGO_10K_raw_qvals.txt", quote = FALSE, sep = "\t")

