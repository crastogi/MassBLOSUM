library(Biostrings)
library(umap)
library(MassBLOSUM)
library(ggplot2)
data("BLOSUM62")
# What is gap cost?
gap_cost = -4

# Define base sequence
base.seq = "HKPQAKSYLPLRLLDY"
seqs = c("HKPQAISYLPYRILDY", "HKPQAKSYLPYRLLDY", "HKPQAKSYLPYRTLDY", "HKPQSKSYLPYRLLDY", "HKPQAKSYLPYRILDY", "YRSPylHHRGGATWQFDY", "DLFRYYYFMWPLDY", "GHYYDIGVFPWDTFDY", "WQQWAGYPRQKYSFDY", "WQQWSGYPRQKYSFDY", "GKSLYGQETTWPHFDY")
# 1 read in sequences, find alignment
pa = pairwiseAlignment(pattern = seqs, subject = base.seq, substitutionMatrix = BLOSUM62)
aligned.seq = as.character(aligned(pa))

# 2 feed this into a Rcpp based vectorizer


# Load 
gifford = read.table("~/Downloads/E1_L14.10000000.tsv.gz", stringsAsFactors = FALSE)
# Set top seq to be sequence with most counts (although should not make much of a difference - maybe some translation/shifting)
base.seq = gifford$V1[which.max(gifford$V5)]
# Create new MassBLOSUM
a = new(MassBLOSUM, BLOSUM62, gap_cost)
a$setBaseSeq(base.seq)

# Now compute 3 sequence sets: R1, R2, R3
R1 = a$computeScores(gifford$V1[gifford$V3>0])
R2 = a$computeScores(gifford$V1[gifford$V4>0])
R3 = a$computeScores(gifford$V1[gifford$V5>0])

# Compute UMAP projection on R1 R2 and R3 (will use them to predict the other datasets)
umap.R1 = umap(R1)
umap.R2 = umap(R2)
umap.R3 = umap(R3)

# Use R2 to predict R3
umap.R3.R2predict = predict(umap.R2, R3)

# Create data frames
df.R2 = data.frame(Seq=gifford$V1[gifford$V4>0], x=umap.R2$layout[,1], y=umap.R2$layout[,2], Counts=gifford$V4[gifford$V4>0], stringsAsFactors = F)
df.R3 = data.frame(Seq=gifford$V1[gifford$V5>0], x=umap.R3$layout[,1], y=umap.R3$layout[,2], Counts=gifford$V5[gifford$V5>0], stringsAsFactors = F)
df.R3.R2 = df.R3
df.R3.R2$x = umap.R3.R2predict[,1]
df.R3.R2$y = umap.R3.R2predict[,2]


ggplot(data=df.R2, aes(x=x, y=y, z=Counts)) +
  stat_summary_hex(fun=sum, binwidth=c(.5,.5)) +
  scale_fill_gradient(expression(paste("Density")), guide=guide_legend(title.hjust = 1), low="#faeffa", high="blue", trans = "log10") +
  theme_bw() +
  theme(aspect.ratio=1, legend.position = c(.95,.05), legend.justification = c("right", "bottom"))

ggplot(data=df.R3, aes(x=x, y=y, z=Counts)) +
  stat_summary_hex(fun=sum, binwidth=c(.5,.5)) +
  scale_fill_gradient(expression(paste("Density")), guide=guide_legend(title.hjust = 1), low="#faeffa", high="blue", trans = "log10") +
  theme_bw() +
  theme(aspect.ratio=1, legend.position = c(.95,.05), legend.justification = c("right", "bottom"))

ggplot(data=df.R2, aes(x=x, y=y, z=Counts)) +
  stat_summary_hex(fun=sum, binwidth=c(.5,.5)) +
  scale_fill_gradient(expression(paste("Density")), guide=guide_legend(title.hjust = 1), low="#faeffa", high="blue", trans = "log10") +
  geom_point(data=df.R3.R2, aes(x=x, y=y), alpha=log10(df.R3.R2$Counts)) + 
  theme_bw() +
  theme(aspect.ratio=1, legend.position = c(.95,.05), legend.justification = c("right", "bottom"))

# Subsample the R1 reads and see what kind of pattern emerges
R1.sub.idx = sample(1:nrow(R1), size = 75000, replace = F)
R1.sub = R1[R1.sub.idx,]
umap.R1.sub = umap(R1.sub)

df.R1 = data.frame(Seq=gifford$V1[gifford$V3>0], Counts=gifford$V3[gifford$V3>0], stringsAsFactors = F)
df.R1 = df.R1[R1.sub.idx,]
df.R1 = cbind(df.R1, x=umap.R1.sub$layout[,1], y=umap.R1.sub$layout[,2])

ggplot(data=df.R1, aes(x=x, y=y, z=Counts)) +
  stat_summary_hex(fun=sum, binwidth=c(.5,.5)) +
  scale_fill_gradient(expression(paste("Density")), guide=guide_legend(title.hjust = 1), low="#faeffa", high="blue", trans = "log10") +
  theme_bw() +
  theme(aspect.ratio=1, legend.position = c(.95,.05), legend.justification = c("right", "bottom"))


giffR3 = gifford$V1[gifford$V5>0]
dm = matrix(data = 0, nrow=length(giffR3), ncol=length(giffR3))
for (i in 1:length(giffR3)) {
  if ((i %% 50)==0) {
    print(i)
  }
  dm[,i] = pairwiseAlignment(pattern = giffR3, subject = giffR3[i], gapOpening=10, gapExtension=4, scoreOnly=T, substitutionMatrix=BLOSUM62)
}
umap.R3.dist = umap(dm, input="dist")
