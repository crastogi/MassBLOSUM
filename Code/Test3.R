c = new(MassBLOSUM, BLOSUM62, 0)
dejan = read.table("~/Desktop/20210222_Anti_streptactin.10000000.tsv.gz", stringsAsFactors = F)

c$setBaseSeq("DDDDDDDDDDD")

R0.dejan = c$computeScores(dejan$V1[dejan$V2>0])
R1.dejan = c$computeScores(dejan$V1[dejan$V3>0])

umap.dejan.R0 = umap(R0.dejan)
umap.dejan.R1 = umap(R1.dejan)

df.R0.dejan = data.frame(Seq=dejan$V1[dejan$V2>0], x=umap.dejan.R0$layout[,1], y=umap.dejan.R0$layout[,2], Counts=dejan$V2[dejan$V2>0], stringsAsFactors = F)
df.R1.dejan = data.frame(Seq=dejan$V1[dejan$V3>0], x=umap.dejan.R1$layout[,1], y=umap.dejan.R1$layout[,2], Counts=dejan$V3[dejan$V3>0], stringsAsFactors = F)

ggplot(data=df.R0.dejan, aes(x=x, y=y, z=Counts)) +
  stat_summary_hex(fun=sum, binwidth=c(.5,.5)) +
  scale_fill_gradient(expression(paste("Density")), guide=guide_legend(title.hjust = 1), low="#faeffa", high="blue", trans = "log10") +
  theme_bw() +
  theme(aspect.ratio=1, legend.position = c(.95,.05), legend.justification = c("right", "bottom"))

ggplot(data=df.R1.dejan, aes(x=x, y=y, z=Counts)) +
  stat_summary_hex(fun=sum, binwidth=c(.5,.5)) +
  scale_fill_gradient(expression(paste("Density")), guide=guide_legend(title.hjust = 1), low="#faeffa", high="blue", trans = "log10") +
  theme_bw() +
  theme(aspect.ratio=1, legend.position = c(.95,.05), legend.justification = c("right", "bottom"))

predict(umap.dejan.R0, c$computeScores("DDDDDDDDDDD"))
