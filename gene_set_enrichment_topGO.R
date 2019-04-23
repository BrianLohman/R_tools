# Load data and setup
library(topGO)
library(org.Hs.eg.db)

# individuals and PRS for stratification
master = read.table("15_Jan_19_Simons_master_ancestry_corrected_PRS.txt", sep = '\t', header = TRUE)

# de novo mutations
denovos = read.table("master_gatk_rufus_med_high_variant_table.txt", sep = '\t', header = TRUE)

# params for analysis
quartile = "first"
trait = "mom_bmi"

# run analysis
if(trait != "mom_bmi"){
      probands = master[master$family_member == 'p1' ,]
  probands = probands[order(probands[paste(trait, "_ancestry_resid", sep = "")]) , ]

    if(quartile == "first"){
            qoi = head(probands, n = 413)
    }

    if(quartile == "fourth"){
            qoi = tail(probands, n = 413)
      }
}

if(trait == "mom_bmi"){
      moms = master[master$family_member == "mo" , ]
  moms = moms[order(moms$NewBMI_ancestry_resid) ,]
    if(quartile == "first"){
            fams = head(moms, n = 413)
      fams = fams$family
          qoi = master[master$family %in% fams & master$family_member == 'p1' ,]

        }
    if(quartile == "fourth"){
            fams = tail(moms, n = 413)
        fams = fams$family
            qoi = master[master$family %in% fams & master$family_member == 'p1' ,]
          }
}

ids = qoi$simons_id

# get variatns of interest and their genes
voi = denovos[denovos$SampleID %in% ids , ]
goi = voi$SYMBOL

# test for enrichment
# prep gene list
x <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
allGenes <- unique(unlist(x))
geneList <- factor(as.integer(allGenes %in% goi))
names(geneList) <- allGenes

# topGO
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize = 1, annot = annFUN.org, mapping = 'org.Hs.eg.db', ID = 'symbol')

# gather results
result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
result_table <- GenTable(GOdata, classic = result, orderBy = "weight", ranksOf = "classic", topNodes = length(score(result)))
result_table$qval = p.adjust(result_table$classic, method = "fdr")

# plot
pdf(file = paste(quartile, "_quartile_topGO_", trait, ".pdf", sep = ""), height = 11, width = 8.5)
par(oma = c(14,1,1,1))
plot(x = seq(1, table(result_table$qval < 0.05)[2], 1),
    y = result_table[result_table$qval < 0.05, "qval"],
    pch =  16, col = "dodgerblue", xaxt = "n",
    xlab = "", ylab = "FDR")
axis(1, at = seq(1, table(result_table$qval < 0.05)[2], 1),
          labels = result_table[result_table$qval < 0.05, "Term"], las = 2)
dev.off()
