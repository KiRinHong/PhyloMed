rm(list = ls())
library(matrixStats)
fileNames <- list.files("./", pattern = ".*.Rdata")
tab <- matrix(NA, nrow = length(fileNames), ncol = 21, dimnames = list(c(),
                                                                       c("CausalType", "N", "A", "B", "Reps",
                                                                         "simes.asym", "fisher.asym", "hmp.asym", 
                                                                         "simes.perm", "fisher.perm", "hmp.perm", 
                                                                         paste0("efdr.",c("js.asym", "jsmix.asym", "js.perm", "jsmix.perm", "sobel"), ".BH"),
                                                                         paste0("dr.",c("js.asym", "jsmix.asym", "js.perm", "jsmix.perm", "sobel"), ".BH"))
))
for (f in 1:length(fileNames)) {
  print(fileNames[f])
  tmp <- as.numeric(sub("[.]$", "", unlist(strsplit(fileNames[f], split = "[aA-zZ]+"))[-1]))
  if(length(tmp) == 4)
    tab[f,1:4] <- tmp
  if(length(tmp) == 5){
    if(tmp[4] < 0){
      tab[f,1:4] <- c(tmp[1:2], tmp[3]*10^tmp[4], tmp[5])
    }
    if(tmp[5] < 0)
      tab[f,1:4] <- c(tmp[1:3], tmp[4]*10^tmp[5])
  }

  load(fileNames[f])
  tab[f, 5] <- nrow(rslt_all$gp)
  tab[f, 6:11] <- round(colMeans(rslt_all$gp[, 1:6] < 0.05, na.rm = T), 4)
  
  tab[f, 12:16] <- round(colMeans(rslt_all$efdr, na.rm = T), 4)
  tab[f, 17:21] <- round(colMeans(rslt_all$dr, na.rm = T), 4)
}

tab <- as.data.frame(tab)
newtab <- tab[order(tab$N, tab$A, tab$B),]
write.csv(newtab, file = "Type1PowerAllCont.csv")
