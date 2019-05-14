# Load data and setup
rm(list = ls())
library(topGO)
library(org.Hs.eg.db)
library(LaCroixColoR)
library(scales)
library(gmodels)

# individuals and PRS for stratification
master = read.table("15_Jan_19_Simons_master_ancestry_corrected_PRS.txt", sep = '\t', header = TRUE)

# de novo mutations
denovos = read.table("master_gatk_rufus_med_high_variant_table.txt", sep = '\t', header = TRUE)

# load simulation data into memory
df = read.table("topGO_10K_raw_qvals.txt", header = TRUE)
colnames(df) = gsub("GO.", "GO:", colnames(df), perl = FALSE)

# function to run analysis and plot
topGO_and_plot = function(quartile, trait, kid){

# get sample IDs of interest
if(trait != "mom_bmi"){
    probands = master[master$family_member == kid ,]
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
          qoi = master[master$family %in% fams & master$family_member == kid ,]
        
    }
    if(quartile == "fourth"){
            fams = tail(moms, n = 413)
        fams = fams$family
            qoi = master[master$family %in% fams & master$family_member == kid ,]
    }
}

ids = qoi$simons_id

# get variatns of interest
voi = denovos[denovos$SampleID %in% ids , ]

goi = voi$SYMBOL

# test for enrichment
x <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")

allGenes <- unique(unlist(x))

geneList <- factor(as.integer(allGenes %in% goi))

names(geneList) <- allGenes

GOdata <- new("topGOdata", 
    ontology = "BP",
    allGenes = geneList,
    nodeSize = 10, 
    annot = annFUN.org,
    mapping = 'org.Hs.eg.db',
    ID = 'symbol'
)

result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

result_table <- GenTable(GOdata, classic = result, orderBy = "weight", ranksOf = "classic", topNodes = length(score(result)))

result_table$qval = p.adjust(result_table$classic, method = "fdr")

write.table(result_table, file = paste(quartile, trait, kid, "data.txt", sep = "_"), sep = "\t", quote = FALSE, row.names = FALSE)
    
# plotting
# collapse to plotting df
plot_df = result_table[result_table$qval <= 0.05 ,]

# check that there are results and quit if not
if(dim(plot_df)[1] == 0){
      stop("no significant results")
}

sim = df[, which(colnames(df) %in% plot_df$GO.ID)]

# convert to log scale
# simulation data
for(i in 1:ncol(sim)){
      sim[,i] = -log(sim[,i])
}

# PRS stratification
plot_df$log_qval = -log(plot_df$qval)

# calculate summary stats
plot_df$sd = NA
plot_df$mean_qval = NA
plot_df$lowerCI = NA
plot_df$upperCI = NA

for(go in plot_df$GO.ID){
    plot_df[plot_df$GO.ID == go, "sd"] = sd(sim[, go])
    plot_df[plot_df$GO.ID == go, "mean_qval"] = mean(sim[, go])
    plot_df[plot_df$GO.ID == go, "lowerCI"] = ci(sim[, go])[2]
    plot_df[plot_df$GO.ID == go, "upperCI"] = ci(sim[, go])[3]
}

# sort based on mean FDR
plot_df = plot_df[order(plot_df$mean_qval) ,]

# plot
pdf(file = paste(quartile,"_quartile_",kid,"_",trait,"_10_topGO_log_observed_vs_10K_simulated.pdf",
    sep = ""), height = 11, width = 8.5
)
par(oma = c(15,1,1,1))
plot(x = seq(1,length(plot_df$GO.ID),1),
    y = plot_df$mean_qval,
    pch = 16,
    col = "dodgerblue",
    xlab = "",
    ylab = "-log q-value (mean and 2 standard deviations)",
    ylim = c(0, max(plot_df$mean_qval) + (max(plot_df$sd) *2)),
    xaxt = "n",
    main = paste(quartile, "quartile", trait, kid, "PRS stratification vs 10K Monte Carlo simulations", sep = " ")
)

# bars for standard deviation
arrows(x0 = seq(1,length(plot_df$GO.ID),1), x1 = seq(1,length(plot_df$GO.ID),1), 
    y0 = plot_df$mean_qval - (plot_df$sd*2), y1 = plot_df$mean_qval + (plot_df$sd*2), 
    angle = 90,
    code = 3,
    col = "dodgerblue",
    length = 0.05
)

# points for PRS stratification
colors = as.numeric(plot_df$log_qval > (plot_df$mean_qval + (plot_df$sd*2))) + 1
points(x = seq(1,length(plot_df$GO.ID),1), y = plot_df$log_qval, pch = 15, col = colors)

axis(1, at = seq(1,length(plot_df$GO.ID),1), 
    labels = paste(paste("(",plot_df$Significant, "/", plot_df$Annotated, ")", sep = ""),
    plot_df$Term, sep = " "), 
    las = 2
)


# add empirical q-value
plot_df$e_pval = NA
for(go in plot_df$GO.ID){
    if(table(df[,go] < plot_df[plot_df$GO.ID == go, "qval"])[1] == 10000){
        plot_df[plot_df$GO.ID == go, "e_pval"] = 0
  } 
  else {
        plot_df[plot_df$GO.ID == go, "e_pval"] = table(df[,go] < plot_df[plot_df$GO.ID == go, "qval"])[2]/10000 
    }
}

text(x = seq(1,length(plot_df$GO.ID),1) + 0.1,
    y = plot_df$mean_qval + (plot_df$sd * 2) + 0.1,
    labels = plot_df$e_pval,
    cex = 0.5,
    srt = 90,
    pos = 3,
    col = as.numeric(plot_df$e_pval < 0.05) + 1
)


legend("topleft", legend = c("Monte Carlo Simulation", "PRS stratification (inside 2 SDs)", "PRS stratification (outside 2 SDs)", "Signficant empirical p-value", "Non-signficant empirical p-value"), 
    col = c("dodgerblue", "black", "red", "red", "black"),
    pch = c(16, 15, 15, 4, 4),
    bty = "n"
)

# plots to explore data
# volcano plot colored py p/q-value
colors = c()
for(i in 1:nrow(result_table)){
      if(as.numeric(result_table[i,6]) > 0.05){
              colors = c(colors, "#C70E7B")
  }
  if(as.numeric(result_table[i,6]) < 0.05){
          if(result_table[i,7] < 0.05){
                    colors = c(colors, "#007BC3")
      } else {
                colors = c(colors, "#009F3F")
          }
    }
}

plot(x = -log(as.numeric(result_table$classic)), 
    y = result_table$Significant/as.numeric(result_table$Expected),
    pch = 16, col = colors,
    xlab =  "-log p-value",
    ylab = "fold enrichment (observed/epected)",
    main = paste(quartile, "quartile", trait, sep = " ")
)

legend("topright", legend = c("p > 0.05", "p < 0.05 but q > 0.05", "q < 0.05"),
    bty = "n",
    col = c("#C70E7B", "#009F3F", "#007BC3"),
    pch = 15
)

# point size is number of genes in pathway
plot(x = -log(as.numeric(result_table$classic)), 
    y = result_table$Significant/as.numeric(result_table$Expected),
    pch = 16, col = alpha("black", 0.1),
    xlab =  "-log p-value",
    ylab = "fold enrichment (observed/epected)",
    main = paste(quartile, "quartile", trait, "\nPoint size = # genes in GO term", sep = " "),
    cex = log(result_table$Annotated)/2
)

# fold enrichment vs number of annotated genes
colors = as.numeric(result_table$qval < 0.05) + 1
pal = c("#FC688200", "#007BC3")
plot(x = result_table$Annotated,
    y = result_table$Significant/as.numeric(result_table$Expected),
    xlab = "Number of genes in GO term",
    ylab = "Fold enrichment (observed/expected)",
    pch = 16,
    col = alpha(pal[colors], 0.5),
    main = paste(quartile, "quartile", trait, sep = " ")
)

legend("topright", legend = c("q > 0.05", "q < 0.05"), col = c("#FC6882", "#007BC3"),
    pch = 15, bty = "n")

# select GO IDs where PRS stratificatino is outside 2 SDs from simulated mean q-value
sig = plot_df[plot_df$log_qval > (plot_df$mean_qval + (plot_df$sd*2)), ]

for(go in sig$GO.ID){
    h = df[, go]
hist(-log(h), col = "grey", breaks = 100, main = sig[sig$GO.ID == go, "Term"], xlab = "-log(Simulated) vs -log(observed) q-values)")
abline(v = plot_df[plot_df$GO.ID == go, "log_qval"], col = "red", lwd = 2)
}

dev.off()
}

# run
for(q in c("first", "fourth")){
  for(t in c("EA", "SCZ", "mom_bmi")){
          for(k in c("p1", "s1")){
                    print(c(q, t, k))
  try(topGO_and_plot(quartile = q, trait = t, kid = k), silent = TRUE)
      }
  }
}
