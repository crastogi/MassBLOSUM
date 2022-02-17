b = new(MassBLOSUM, BLOSUM62, 0)
all.seqs = c("HKPQAISYLPYRILDY", "HKPQAKSYLPYRLLDY", "HKPQAKSYLPYRTLDY", "HKPQSKSYLPYRLLDY", "HKPQAKSYLPYRILDY", "YRSPHHRGGATWQFDY", "GHYYDIGVFPWDTFDY", "WQQWAGYPRQKYSFDY", "WQQWSGYPRQKYSFDY", "GKSLYGQETTWPHFDY")

# 1 select a base sequence
b$setBaseSeq(all.seqs[1])
# 2 compute vector values
pos1 = b$computeScores(all.seqs)
dist(pos1, method = "euclidean", upper = T)

# 1 select a base sequence
b$setBaseSeq(all.seqs[5])
pos2 = b$computeScores(all.seqs)
dist(pos2, method = "euclidean", upper = T)


dist = function(inMatrix) {
  for (i in 1:nrow)
}


M = BLOSUM62
M = -M
outM = M
for (i in 1:nrow(M)) {
  for (j in 1:nrow(M)) {
    outM[i,j] = M[i,j] - min(M[i,i], M[j,j])
  }
}

for (i in 1:nrow(M)) {
  for (j in 1:nrow(M)) {
    dxz = outM[i,j]
    for (k in 1:nrow(M)) {
      if (dxz > outM[i,k]+outM[k,j]) {
        print(paste(i,j,k))
      }
    }
  }
}

for (i in 1:nrow(M)) {
  for (j in 1:nrow(M)) {
    dxz = outM[i,j]
    for (k in 1:nrow(M)) {
      if (dxz != outM[i,k]+outM[k,j]) {
        print(paste(rownames(outM)[i],rownames(outM)[j],rownames(outM)[k]))
      }
    }
  }
}

for (i in 1:nrow(M)) {
  for (j in 1:nrow(M)) {
    dxz = M[i,j]
    for (k in 1:nrow(M)) {
      if (dxz > M[i,k]+M[k,j]) {
        print(paste(i,j,k))
      }
    }
  }
}