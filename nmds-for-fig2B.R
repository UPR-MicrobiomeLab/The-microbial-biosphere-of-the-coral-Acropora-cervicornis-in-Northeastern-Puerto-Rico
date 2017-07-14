library(vegan)

data= read.csv("38distmatrix.csv", header=TRUE, row.names=1)

meta = read.table("metacoral.txt", sep="\t", row.names=1, header=T)

metaMDS_coral = metaMDS(t(data), autotransform = F)

pdf(file ="nmds_coral-depth.pdf")

plot(metaMDS_coral$points, pch=meta$coral, col=meta$depth,lwd=5, main="NMDS by coral depth”, cex.lab = 1.2, cex.axis = 1.2, cex.main =1.2 ,xlab="Axis 1", ylab="Axis 2")

text(metaMDS_coral$points, labels = rownames(meta), pos=1, cex=0.5)

metaMDS_coral$stress

dev.off()

