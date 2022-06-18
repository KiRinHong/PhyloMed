rm(list = ls())
library(matrixStats)
setwd("~/Documents/Project/PhyloMed/Simulation/PropOfNull/Result")
fileNames <- list.files("./", pattern = ".*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 10, dimnames = list(c(),
                                                                       c("CausalType", "N", "A", "B",
                                                                         "Pi00.prod-true", "Pi00.minus-true",
                                                                         "Pi10.prod-true", "Pi10.minus-true",
                                                                         "Pi01.prod-true", "Pi01.minus-true")))

for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tab[f, 1:4] <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  load(fileNames[f])
  tmp <- sweep(rslt$alpha00.mat[,-1], 1, rslt$alpha00.mat[,1], "-")
  tab[f, 5:6] <- paste0(round(colMeans(tmp), 3), " (", round(colSds(tmp), 3), ") ")
  tmp <- sweep(rslt$alpha10.mat[,-1], 1, rslt$alpha10.mat[,1], "-")
  tab[f, 7:8] <- paste0(round(colMeans(tmp), 3), " (", round(colSds(tmp), 3), ") ")
  tmp <- sweep(rslt$alpha01.mat[,-1], 1, rslt$alpha01.mat[,1], "-")
  tab[f, 9:10] <- paste0(round(colMeans(tmp), 3), " (", round(colSds(tmp), 3), ") ")
}

tab <- as.data.frame(tab)
# large effect size A
tab_sel_l <- tab[c(3,4,5,8,9),]
# small effect size A
tab_sel_s <-  tab[c(1,2,6,7),]
write.csv(tab, file = "PropNull.csv")
