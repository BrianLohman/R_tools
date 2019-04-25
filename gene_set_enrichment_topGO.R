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


# when running in conjuctin with simulation of null distribution
# collapse to plotting df
plot_df = result_table[result_table$qval <= 0.05, ]
colnames(df)[1] = "Term"
sim = df[df$Term %in% plot_df$Term , ]

plot_df = merge(plot_df, sim, by = "Term", all.x = TRUE)

# add +/- 1 SD
plot_df$upper = plot_df$meanFDR + plot_df$sdFDR
plot_df$lower = plot_df$meanFDR - plot_df$sdFDR

# sort based on mean FDR
plot_df = plot_df[order(plot_df$meanFDR) ,]

# plot
pdf(file = paste(quartile, "_quartile", trait, "_topGO_observed_vs_simulated.pdf", sep = ""), height = 11, width = 8.5)
par(oma = c(14,1,1,1))
plot(x = seq(1,length(plot_df$Term),1),
     y = plot_df$meanFDR,
     pch = 16,
     col = "dodgerblue",
     xlab = "",
     ylab = "FDR (mean and standard deviation)",
     ylim = c(0, max(plot_df$meanFDR) + max(plot_df$sdFDR)),
     xaxt = "n",
     main = paste(quartile, "quartile", trait, "PRS stratification vs Monte Carlo simulation", sep = " ")
)

axis(1, at = seq(1,length(plot_df$Term),1), labels = plot_df$Term, las = 2)

arrows(x0 = seq(1,length(plot_df$Term),1), x1 = seq(1,length(plot_df$Term),1), 
       y0 = plot_df$lower, y1 = plot_df$upper, 
       angle = 90,
       code = 3,
       col = "dodgerblue",
       length = 0.1
)

points(x = seq(1,length(plot_df$Term),1), y = plot_df$qval, pch = 15, col = "firebrick")

legend("topleft", legend = c("Monte Carlo Simulation", "PRS stratification"), 
       col = c("dodgerblue", "firebrick"),
       pch = c(16, 15),
       bty = "n"
)

dev.off()
