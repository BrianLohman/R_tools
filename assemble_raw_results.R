rm(list = ls())

# glob in RData files
filenames <- Sys.glob("*.RData")

# empty sapce to store master resutls list and iteration ID
master_results = list()
iterations = c()
sample_ids = list()
genes = list()

# unpack RData objects
for(i in filenames){
    # load file into  memory
    load(i)

    # append to master list
    master_results = c(master_results, mc_results)

    # add file name to list of iterations
    iteration = c(iterations, strsplit(filenames[1], split = "_")[[1]][1])

    # add samples in the iteration to the master sample list
    sample_ids = c(sample_ids, ids)
    
    # get genes for each iteration and save to master
    genes = c(genes, as.character(goi))
}

# check length 10000
if(length(master_results) != 2){
    stop("incorrect number of simulations loaded")
}

# get list of GO terms
uniq_go = master_results[[1]][,1]

# print number of go terms included in analysis
length(uniq_go)

# make table of iteration and sample IDs
sample_dict = as.data.frame(matrix(nrow = 413, ncol = 2))
colnames(sample_dict) = iterations
for(i in 1:2){
    sample_dict[,i] = sample_ids[[i]][1]
}
sample_dict

# save to file
write.table(sample_dict, file = "topGO_10K_sample_dict.txt", quote = FALSE, sep = "\t")

# make table of genes for each iteration
gene_table = as.data.frame(matrix(nrow = 2, ncol = length(genes[[1]])))
row.names(gene_table) = iterations
colnames(gene_table) = genes[[1]]

for(i in row.names(gene_table)){
    

}

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

