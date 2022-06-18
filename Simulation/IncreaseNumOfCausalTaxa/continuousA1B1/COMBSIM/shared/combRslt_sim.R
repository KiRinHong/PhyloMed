rm(list=ls())
args = (commandArgs(trailingOnly=TRUE))
if(length(args) == 4){
  causal.type = as.numeric(args[1]) # can take values 1,2,3 corresponding to four approaches to simulate causal species
  n.sample = as.numeric(args[2]) # number of samples
  A = as.numeric(args[3]) # correspond to A in the document, control the perturbation level for species S_ab and S_a in treatment group
  B = as.numeric(args[4]) # correspond to B in the document, control the mediation effect in the outcome model for species S_ab and S_b
} else {
  cat('usage: Rscript combRslt_sim.R <causalType> <n.sample> <A> <B>\n', file = stderr())
  stop()
}
fileNames <- list.files(pattern = paste0("^compType", causal.type, "Nsample", n.sample, ".*", "A", A, "B", B, ".Rdata$"))
rslt_all <- list()
for (f in 1:length(fileNames)) {
  print(fileNames[f])
  load(fileNames[f])
  if(length(rslt_all) == 0){
    rslt_all = rslt
  }else{
    for (i in 1:length(rslt_all)) {
      if(is.null(rslt_all[[i]])){
        next
      }else{
        rslt_all[[i]] = rbind(rslt_all[[i]], rslt[[i]])
      }
    }
  }
}
filename = paste0("allType", causal.type, "Nsample", n.sample, "A", A, "B", B, ".Rdata")
save(rslt_all, file = filename)



