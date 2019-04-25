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
x <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
allGenes <- unique(unlist(x))

# setup space for results
mc_results = list()

# run
for(i in seq(1,5000,1)){
    # get sample ids (random)
    poi = sample(proband_ids, 413)
    qoi = probands[probands$simons_id %in% poi , ]  
    ids = qoi$simons_id

    # get genes which harbor variants of interest
    voi = denovos[denovos$SampleID %in% ids , ]
    goi = voi$SYMBOL
    geneList <- factor(as.integer(allGenes %in% goi))
    names(geneList) <- allGenes

    # run enrichment test
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize = 1, annot = annFUN.org, mapping = 'org.Hs.eg.db', ID = 'symbol')
    result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

    # gather results
    result_table <- GenTable(GOdata, classic = result, orderBy = "weight", ranksOf = "classic", topNodes = length(score(result)))

    # FDR correction
    result_table$qval = p.adjust(result_table$classic, method = "fdr")

    mc_results[[i]] = result_table

}

# save
save(mc_results, file = "gene_set_enrichment_Monte_Carlo_results.RData")

# parse results and plot
load("gene_set_enrichment_Monte_Carlo_results.RData")

# make list of uniq go terms
go_terms = c()
for(i in 1:length(mc_results)){
      go_terms = c(go_terms, mc_results[[i]][,2])
}
uniq_go = unique(go_terms)

# loop over go terms, calculating mean and sd
df = as.data.frame(matrix(nrow = length(uniq_go), ncol = 3))
colnames(df) = c("go_term", "meanFDR", "sdFDR")
df$go_term = uniq_go

for(go in uniq_go){
      fdr = c()
  for(i in 1:length(mc_results)){
          tmp = mc_results[[i]]
      fdr = c(fdr, tmp[tmp$Term == go, "qval"])
        }
    df[df$go_term == go, "meanFDR"] = mean(fdr)
    df[df$go_term == go, "sdFDR"] = sd(fdr)
}

pdf(file = "test_simulation_plot.pdf", height = 11, width = 8.5)
par(oma = c(14,1,1,1))
plot(x = seq(1, table(df$meanFDR < 0.05)[2], 1), 
          y = df[df$meanFDR < 0.05, "meanFDR"],
               pch =  16, col = "dodgerblue", xaxt = "n",
                    xlab = "", ylab = "FDR (mean and SD)")
axis(1, at = seq(1, table(df$meanFDR < 0.05)[2], 1), 
          labels = df[df$meanFDR < 0.05, "go_term"], las = 2)
dev.off()
