rm(list = ls())
library(matrixStats)
fileNames <- list.files("./", pattern = ".*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 26, dimnames = list(c(),
                                                                       c("CausalType", "N", "A", "B", "Reps",
                                                                         "dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.clr.MedTest", "dist.omn.MedTest",
                                                                         "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA", "dist.clr.MODIMA", "dist.bonf.MODIMA",
                                                                         "ldm.med.global",
                                                                         "simes.asym", "fisher.asym", "hmp.asym", 
                                                                         "simes.perm", "fisher.perm", "hmp.perm", 
                                                                         "ccmm.boot.tide", "ccmm.norm.tide")
))
for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tab[f, 1:4] <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  load(fileNames[f])
  tab[f, 5] <- nrow(rslt_all$gp)
  tab[f, 6:16] <- round(colMeans(rslt_all$gp[, 1:11] < 0.05, na.rm = T), 4)
  tab[f, 17] <- round(mean(rowMins(cbind(rslt_all$gp[, 7:11] * 4, 1), na.rm = T)  < 0.05), 4) 
  tab[f, 18:26] <- round(colMeans(rslt_all$gp[, 12:20] < 0.05, na.rm = T), 4)
}

tab <- as.data.frame(tab)
newtab <- tab[order(tab$N, tab$A, tab$B),]
write.csv(newtab, file = "PowerAllCont.csv")
