rm(list = ls())
library(matrixStats)
fileNames <- list.files("./", pattern = ".*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 42, dimnames = list(c(),
                                                                       c("CausalType", "N", "A", "B",
                                                                         "dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.omn.MedTest",
                                                                         "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA", "dist.bonf.MODIMA",
                                                                         "simes.asym", "fisher.asym", "hmp.asym", 
                                                                         "simes.perm", "fisher.perm", "hmp.perm", 
                                                                         "ccmm.boot.tide", "ccmm.norm.tide",
                                                                         paste0("efdr.",c("js.asym", "jsmix.asym", "js.perm", "jsmix.perm", "sobel"), ".BH"),
                                                                         paste0("efdr.",c("js.asym", "jsmix.asym", "js.perm", "jsmix.perm", "sobel"), ".HMP"),
                                                                         paste0("dr.",c("js.asym", "jsmix.asym", "js.perm", "jsmix.perm", "sobel"), ".BH"),
                                                                         paste0("dr.",c("js.asym", "jsmix.asym", "js.perm", "jsmix.perm", "sobel"), ".HMP"))
))
for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tab[f, 1:4] <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  load(fileNames[f])
  tab[f, 5:13] <- round(colMeans(rslt_all$gp[, 1:9] < 0.05), 4)
  tab[f, 14] <- round(mean(rowMins(cbind(rslt_all$gp[, 6:9] * 4, 1))  < 0.05), 4) 
  tab[f, 15:22] <- round(colMeans(rslt_all$gp[, 10:17] < 0.05), 4)
  
  if(all(tab[f, 3:4] == 1)){
    tab[f, 23:32] <- round(colMeans(rslt_all$efdr), 4)
    tab[f, 33:42] <- round(colMeans(rslt_all$dr), 4)
  }
  
}

tab <- as.data.frame(tab)
newtab <- tab[order(tab$N, tab$A, tab$B),]
write.csv(newtab, file = "Type1PowerAllCont.csv")
