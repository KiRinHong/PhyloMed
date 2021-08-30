rm(list = ls())
# partial overlap
tmp <- do.call(expand.grid, list(type=1:2, nsample=c(50,200), sim=1:2000, A=1, B=1))
tmp$outfile <- paste0("compType", tmp$type, "Nsample", tmp$nsample, 
                      "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

tmp <- do.call(expand.grid, list(type=1:2, nsample=c(50,200), A=1, B=1))
tmp$outfile <- paste0("allType", tmp$type, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

# small effect size
tmp <- do.call(expand.grid, list(type=1:2, nsample=200, sim=1:2000, A=c(0,0.001), B=c(0,0.5)))
tmp <- tmp[-which(tmp$A == 0 & tmp$B == 0),]
tmp$outfile <- paste0("compType", tmp$type, "Nsample", tmp$nsample, 
                      "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

tmp <- do.call(expand.grid, list(type=1:2, nsample=200, A=c(0,0.001), B=c(0,0.5)))
tmp <- tmp[-which(tmp$A == 0 & tmp$B == 0),]
tmp$outfile <- paste0("allType", tmp$type, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)


# increase num of mediating taxa
tmp <- do.call(expand.grid, list(type=c(1,2,5), nsample=c(50,200), sim=1:2000, A=0.001, B=0.5))
tmp$outfile <- paste0("compType", tmp$type, "Nsample", tmp$nsample, 
                      "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

tmp <- do.call(expand.grid, list(type=c(1,2,5), nsample=c(50,200), A=0.001, B=0.5))
tmp$outfile <- paste0("allType", tmp$type, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
