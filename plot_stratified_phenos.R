prs_of_interest = c("Educational.Attainment_ancestry_resid_19", "Autism.Spectrum.Disorder_ancestry_resid_19", "ADHD.2017_ancestry_resid_19", "Educational.Attainment_ancestry_resid_19", "Major.Depressive.Disorder_ancestry_resid_19", "Schizophrenia_ancestry_resid_19", "Body.Mass.Index_ancestry_resid_19")

for(p in prs_of_interest){
      prs = p
  
  if(prs == "Body.Mass.Index_ancestry_resid_19"){
      moms = master[master$family_member == "mo" , ]
      moms = moms[order(moms$Body.Mass.Index_ancestry_resid_19) , ]
      fourth_fams = tail(moms$family, n = 325)
      first_fams = head(moms$family, n = 325)
              
      lower = master[master$family_member == "p1" & master$family %in% first_fams , ]
      upper = master[master$family_member == "p1" & master$family %in% fourth_fams , ]
               
      } else {
      
      probands = master[master$family_member == "p1" ,]
      sorted = probands[order(probands[. which(colnames(probands) == p)]), ]
      fourth = tail(sorted$IID, n = 325)
      first = head(sorted$IID, n = 325)
      lower = master[master$IID %in% first ,]
      upper = master[master$IID %in% fourth ,]
                        
    }
    traits = c("bmi", "head_circumference", "srs_total", "scq_life_total", "father_age_birth", "sex","gca_standard" )
    
    pdf(file = paste(prs, ".pdf", sep = ""), height = 8, width = 11)
    par(mfrow = c(3,3))
    for(t in traits) {
        y = lower[, t]
        z = upper[, t]
              
        # fathers age at bith is 0 for many probands (??)
        if(t == "father_age_birth"){
            y = y[y > 0]
            z = z[z > 0]
        }
              
        if(t == "sex"){
            stat1 = as.numeric(binom.test(table(y), p = 0.87)[3])
            stat2 = as.numeric(binom.test(table(z), p = 0.87)[3])
            stat = paste(stat1, stat2, sep = " p = ")
                        
            } else {
                stat = t.test(y, z, var.equal = FALSE)[3]
            }
                  
            plot(x = c(jitter(rep(1.2,length(y))), jitter(rep(1.8, length(z)))), 
                y = c(y,z), pch = 16, col = c(rep("dodgerblue", length(y)), rep("firebrick", length(z))), 
                xlim = c(1,2), 
                xaxt = "n", 
                xlab = paste( "Quartile", paste("p = ", stat, sep = ""), sep = "\n"),
                ylab = paste(t, "(mean +/- 1 SD)", sep = " "),
                main = paste(sub("_ancestry_resid_19", "", prs), "PRS", sep = " "),
                                                                                 )
                axis(side = 1, at = c(1.2, 1.8), labels = c("First", "Fourth"))
                points(x = 1.3, y = mean(y, na.rm = TRUE), col = "dodgerblue", pch = 15)
                points(x = 1.9, y = mean(z, na.rm = TRUE), col = "firebrick", pch = 15)
                arrows(x0 = 1.3, y0 = mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), x1 = 1.3 ,y1 = mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), 
                    angle = 90, code = 3, length = 0.1)
                arrows(x0 = 1.9, y0 = mean(z, na.rm = TRUE) - sd(z, na.rm = TRUE), x1 = 1.9 ,y1 = mean(z, na.rm = TRUE) + sd(z, na.rm = TRUE), 
                    angle = 90, code = 3, length = 0.1)
    }
    dev.off()
}
