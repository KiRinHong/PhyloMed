rm(list = ls())
fileNames <- list.files("./", pattern = ".*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 31, dimnames = list(c(),
                                                                       c("CausalType", "N", "A", "B",
                                                                         "dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.omn.MedTest",
                                                                         "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA",
                                                                         "simes.asym", "fisher.asym", "hmp.asym", 
                                                                         "simes.perm", "fisher.perm", "hmp.perm", 
                                                                         "ccmm.boot.tide", "ccmm.norm.tide",
                                                                         paste0("efdr.",c("js.asym", "jsmix.asym", "js.perm", "jsmix.perm", "sobel")),
                                                                         paste0("dr.",c("js.asym", "jsmix.asym", "js.perm", "jsmix.perm", "sobel")))
))
for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tab[f, 1:4] <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  load(fileNames[f])
  tab[f, 5:21] <- round(colMeans(rslt_all$gp < 0.05, na.rm = T), 4)
  if(all(tab[f, 3:4] == 1)){
    tab[f, 22:26] <- round(colMeans(rslt_all$efdr), 4)
    tab[f, 27:31] <- round(colMeans(rslt_all$dr), 4)
  }
  
}

tab <- as.data.frame(tab)
newtab <- tab[order(tab$N, tab$A, tab$B),]
write.csv(newtab, file = "Type1PowerAllCont.csv")
