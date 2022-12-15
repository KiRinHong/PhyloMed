rm(list = ls())
setwd("~/Documents/Project/PhyloMed/Simulation/")
# regular 
tmp <- do.call(expand.grid, list(type=1:2, nsample=c(50,200), sim=1:2000, A=c(0,0.5), B=c(0,0.5)))
tmp <- tmp[-which(tmp$type == 2 & tmp$A == 0 & tmp$B == 0),]
tmp$outfile <- paste0("compType", tmp$type, "Nsample", tmp$nsample, 
                      "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "ProposedModel_Main/continuous_JC/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
write.table(tmp, file = "ProposedModel_Main/binary_JC/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

tmp <- do.call(expand.grid, list(type=1:2, nsample=c(50,200), A=c(0,0.5), B=c(0,0.5)))
tmp <- tmp[-which(tmp$type == 2 & tmp$A == 0 & tmp$B == 0),]
tmp$outfile <- paste0("allType", tmp$type, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "ProposedModel_Main/continuous_JC/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
write.table(tmp, file = "ProposedModel_Main/binary_JC/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

# # small effect size
# tmp <- do.call(expand.grid, list(type=1:2, nsample=200, sim=1:2000, A=c(0,0.001), B=c(0,0.5)))
# tmp <- tmp[-which(tmp$A == 0 & tmp$B == 0),]
# tmp$outfile <- paste0("compType", tmp$type, "Nsample", tmp$nsample, 
#                       "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
# anyDuplicated(tmp)
# write.table(tmp, file = "SmallEffectSize/continuous/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
# 
# tmp <- do.call(expand.grid, list(type=1:2, nsample=200, A=c(0,0.001), B=c(0,0.5)))
# tmp <- tmp[-which(tmp$A == 0 & tmp$B == 0),]
# tmp$outfile <- paste0("allType", tmp$type, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
# anyDuplicated(tmp)
# write.table(tmp, file = "SmallEffectSize/continuous/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
# 
# # partial overlap
# tmp <- do.call(expand.grid, list(type=1:2, nsample=c(50,200), sim=1:2000, A=1, B=1))
# tmp$outfile <- paste0("compType", tmp$type, "Nsample", tmp$nsample, 
#                       "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
# anyDuplicated(tmp)
# write.table(tmp, file = "PartialOverlap/continuous/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
# write.table(tmp, file = "PartialOverlap/binary/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
# 
# tmp <- do.call(expand.grid, list(type=1:2, nsample=c(50,200), A=1, B=1))
# tmp$outfile <- paste0("allType", tmp$type, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
# anyDuplicated(tmp)
# write.table(tmp, file = "PartialOverlap/continuous/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
# write.table(tmp, file = "PartialOverlap/binary/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

# randomly scattered, increase num of mediating taxa
tmp <- do.call(expand.grid, list(type=c(1,2,5), nsample=c(50,200), sim=1:2000, A=c(0.001,0.01,0.1,0.5), B=0.5))
tmp$outfile <- paste0("compType", tmp$type, "Nsample", tmp$nsample, 
                      "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "IncreaseNumOfCausalTaxa/continuous/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

tmp <- do.call(expand.grid, list(type=c(1,2,5), nsample=c(50,200), A=c(0.001,0.01,0.1,0.5), B=0.5))
tmp$outfile <- paste0("allType", tmp$type, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "IncreaseNumOfCausalTaxa/continuous/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

# rare taxa 
tmp <- do.call(expand.grid, list(type=c(1,2), nsample=c(50,200), sim=1:2000, A=0.5, B=0.5))
tmp$outfile <- paste0("compType", tmp$type, "Nsample", tmp$nsample, 
                      "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "RareTaxa/continuous/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
write.table(tmp, file = "RareTaxa/binary/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

tmp <- do.call(expand.grid, list(type=c(1,2), nsample=c(50,200), A=0.5, B=0.5))
tmp$outfile <- paste0("allType", tmp$type, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "RareTaxa/continuous/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
write.table(tmp, file = "RareTaxa/binary/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

# Increase effect size --> empirical FDR
tmp <- do.call(expand.grid, list(type=1:3, nsample=200, sim=1:2000, A=c(0,0.1,0.01,1e-3,1e-4,1e-5), B = 0.5))
tmp$outfile <- paste0("compType", tmp$type, "Nsample", tmp$nsample, 
                      "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "IncreaseEffectSize/continuous/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

tmp <- do.call(expand.grid, list(type=1:3, nsample=200, A=c(0,0.1,0.01,1e-3,1e-4,1e-5), B = 0.5))
tmp$outfile <- paste0("allType", tmp$type, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "IncreaseEffectSize/continuous/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

# pseduo count
tmp <- do.call(expand.grid, list(type=1:2, count=1:2, nsample=c(50,200), sim=1:2000, A=c(0,0.5), B=c(0,0.5)))
tmp <- tmp[-which(tmp$type == 2 & tmp$A == 0 & tmp$B == 0),]
tmp$outfile <- paste0("compType", tmp$type, "Count", tmp$count, "Nsample", tmp$nsample, 
                      "Sim", tmp$sim, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "PseudoCount/continuous/GETSIM/input.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
tmp <- do.call(expand.grid, list(type=1:2, count=1:2, nsample=c(50,200), A=c(0,0.5), B=c(0,0.5)))
tmp <- tmp[-which(tmp$type == 2 & tmp$A == 0 & tmp$B == 0),]
tmp$outfile <- paste0("allType", tmp$type, "Count", tmp$count, "Nsample", tmp$nsample, "A", tmp$A, "B", tmp$B, ".Rdata")
anyDuplicated(tmp)
write.table(tmp, file = "PseudoCount/continuous/COMBSIM/comb_sim.txt", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

