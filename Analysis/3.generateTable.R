# significance level 0.01

rm(list = ls())
setwd("~/Project/QH0020/medtest/PhyloMed/Analysis/")
### continous outcome
fileNames <- list.files("../Simulation/continuous_JC/OUTPUT/shared/Type/", pattern = ".*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 21, dimnames = list(c(),
                                                                       c("CausalType", "N", "A", "B",
                                                                         "dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.omn.MedTest",
                                                                         "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA",
                                                                         "simes.asym", "fisher.asym", "hmp.asym", 
                                                                         "simes.perm", "fisher.perm", "hmp.perm", 
                                                                         "ccmm.boot.tide", "ccmm.norm.tide")
))
for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tmp <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  if(tmp[3] == 1 & tmp[4] == 1) next
  tab[f, 1:4] <- tmp
  load(paste0("../Simulation/continuous_JC/OUTPUT/shared/Type/", fileNames[f]))
  tab[f, 5:21] <- round(colMeans(rslt_all$gp < 0.01, na.rm = T), 4)
}

tab <- as.data.frame(tab[complete.cases(tab[,1:19]),])
write.csv(tab, file = "../Simulation/continuous_JC/RESULT/Type1PowerAllCont_01.csv")

### binary outcome

fileNames <- list.files("../Simulation/binary_JC/OUTPUT/shared/Type/", pattern = ".*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 21, dimnames = list(c(),
                                                                       c("CausalType", "N", "A", "B",
                                                                         "dist.bc.MedTest", "dist.jac.MedTest", "dist.uw.MedTest", "dist.w.MedTest", "dist.omn.MedTest",
                                                                         "dist.bc.MODIMA", "dist.jac.MODIMA", "dist.uw.MODIMA", "dist.w.MODIMA",
                                                                         "simes.asym", "fisher.asym", "hmp.asym", 
                                                                         "simes.perm", "fisher.perm", "hmp.perm", 
                                                                         "ccmm.boot.tide", "ccmm.norm.tide")
))
for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tmp <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  if(tmp[3] == 1 & tmp[4] == 1) next
  tab[f, 1:4] <- tmp
  load(paste0("../Simulation/binary_JC/OUTPUT/shared/Type/", fileNames[f]))
  tab[f, 5:21] <- round(colMeans(rslt_all$gp < 0.01, na.rm = T), 4)
}

tab <- as.data.frame(tab[complete.cases(tab[,1:19]),])
write.csv(tab, file = "../Simulation/binary_JC/RESULT/Type1PowerAllCont_01.csv")
