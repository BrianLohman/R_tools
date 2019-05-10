rm(list = ls())

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
        master_samples[[100*k+j]] = as.character(ids[[j]])
    
        # get genes for each iteration and save to master
        master_genes[[100*k+j]] = as.character(goi[[j]])
    }
}

# check length 10000
#if(length(master_results) != 10000){
#    stop("incorrect number of simulations loaded")
#}

# get list of GO terms
uniq_go = master_results[[1]][,1]

# print number of go terms included in analysis
length(uniq_go)

# make table of iteration and sample IDs
sample_dict = as.data.frame(matrix(nrow = 413, ncol = 200))
colnames(sample_dict) = iterations
for(i in 1:2){
    sample_dict[,i] = master_samples[[i]][1]
}
sample_dict

# save to file
write.table(sample_dict, file = "topGO_10K_sample_dict.txt", quote = FALSE, sep = "\t")

# make table of genes for each iteration
# maybe replace this with a different approach:
#   convert each element in the list to a data frame with 1 row
#   make the column names the gene names
#   merge onto prior element in the list

gene_table = as.data.frame(matrix(nrow = 2, ncol = length(master_genes[[1]])))
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

