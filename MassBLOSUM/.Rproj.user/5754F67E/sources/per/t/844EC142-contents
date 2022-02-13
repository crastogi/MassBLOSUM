# Define base sequence and blosum matrix
base.seq = "HKPQAKSYLPLRLLDY"
bm = BLOSUM62
gapCost = -15

a = new(MassBLOSUM, bm, gapCost)
a$setBaseSeq(base.seq)

# Validate matrix
alpha = colnames(bm)

output = NULL
for (curr.char in unlist(strsplit(base.seq, split = ""))) {
  output = cbind(output, bm[,which(curr.char==alpha)])
}
output = rbind(output,gapCost)


# Define seqs
seqs = c("HKPQAISYLPYRILDY", "HKPQAKSYLPYRLLDY", "HKPQAKSYLPYRTLDY", "HKPQSKSYLPYRLLDY", "HKPQAKSYLPYRILDY", "YRSPHHRGGATWQFDY", "DLFRYYYFMWPLDY", "GHYYDIGVFPWDTFDY", "WQQWAGYPRQKYSFDY", "WQQWSGYPRQKYSFDY", "GKSLYGQETTWPHFDY")
a$computeScores(seqs)
