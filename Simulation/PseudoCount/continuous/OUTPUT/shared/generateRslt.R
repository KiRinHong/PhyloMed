rm(list = ls())
library(matrixStats)
fileNames <- list.files("./", pattern = ".*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 12, dimnames = list(c(),
                                                                       c("CausalType", "CountType", "N", "A", "B", "Reps",
                                                                         "simes.asym", "fisher.asym", "hmp.asym", 
                                                                         "simes.perm", "fisher.perm", "hmp.perm")
))
for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tab[f, 1:5] <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  load(fileNames[f])
  tab[f, 6] <- nrow(rslt_all$gp)
  tab[f, 7:12] <- round(colMeans(rslt_all$gp[, 1:6] < 0.05, na.rm = T), 4)
}

tab <- as.data.frame(tab)
newtab <- tab[order(tab$N, tab$A, tab$B),]
write.csv(newtab, file = "Type1PowerAllCont.csv")
