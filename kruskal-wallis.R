args <- commandArgs(trailingOnly = TRUE)
options(warn = -1)

expr = read.csv(args[1],header=FALSE)
A = as.numeric(expr[1,])
B = as.numeric(expr[2,])
S = (A/(A+B))
S[is.nan(S)] = 0.5

cell_cycle_stage = read.csv(args[2],header=FALSE,row.names=1)
G1_inds = which(cell_cycle_stage[1] == "G1")
S_inds = which(cell_cycle_stage[1] == "S")

G2M_inds = which(cell_cycle_stage[1] == "G2M")
cat(kruskal.test(list(S[G1_inds],S[S_inds],S[G2M_inds]))$p.value)