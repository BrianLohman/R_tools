# install necessary packages
#if (!requireNamespace("BiocManager"))
#        install.packages("BiocManager")
#BiocManager::install()

#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install("topGO", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db", version = "3.8")

library(topGO)
library(org.Hs.eg.db)

# individuals and PRS for stratification
master = read.table("15_Jan_19_Simons_master_ancestry_corrected_PRS.txt", sep = '\t', header = TRUE)
probands = master[master$family_member == 'p1' ,]
proband_ids = probands$simons_id

# de novo mutations
denovos = read.table("master_gatk_rufus_med_high_variant_table.txt", sep = '\t', header = TRUE)

# gene symbol to GO term mapping
x = annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
allGenes = unique(unlist(x))

# setup space for results and to track samples and genes in each iteration
mc_results = list()
genes = list()
samples = list()

# set seed (fractions of second at current time * PID)
s = as.numeric(format(Sys.time(), "%OS3")) * 1000
s
pid = Sys.getpid()
pid
set.seed(s*as.numeric(pid))

# run
for(i in seq(1,100,1)){
    # get sample ids (random)
    poi = sample(proband_ids, 413)
    qoi = probands[probands$simons_id %in% poi , ]  
    ids = qoi$simons_id

    # get genes which harbor variants of interest
    voi = denovos[denovos$SampleID %in% ids , ]
    goi = voi$SYMBOL
    geneList = factor(as.integer(allGenes %in% goi))
    names(geneList) = allGenes

    # run enrichment test
    GOdata = new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize = 10, annot = annFUN.org, mapping = 'org.Hs.eg.db', ID = 'symbol')
    result = runTest(GOdata, algorithm = "classic", statistic = "fisher")

    # gather results
    result_table = GenTable(GOdata, classic = result, orderBy = "weight", ranksOf = "classic", topNodes = length(score(result)))

    # FDR correction
    result_table$qval = p.adjust(result_table$classic, method = "fdr")

    mc_results[[i]] = result_table
    genes[[i]] = as.character(goi)
    samples[[i]] = as.character(ids)
}

# save file as random seed * finish time
s2 = as.numeric(format(Sys.time(), "%OS3")) * 1000
save(mc_results, samples, genes, file = paste(as.character(s*s2), "_gene_set_enrichment_Monte_Carlo_results.RData", sep = ""))
