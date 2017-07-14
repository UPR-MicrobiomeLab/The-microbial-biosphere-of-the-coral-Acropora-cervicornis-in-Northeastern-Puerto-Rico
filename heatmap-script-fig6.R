 library(vegan)

library(devtools)

meta <- read.csv("Coral_Mapping.csv",header=TRUE,row.names=1, sep=",")
 
dat1 <- read.csv("gtest38.csv", header=TRUE,row.names=1, sep=",")
 
dat3<-t(dat1)
 
data.prop_relA<- data.matrix(dat3)

data.prop_relA<- data.matrix(dat3)
 
mydist=function(c) {vegdist(c,method="euclidean")}
 
myclust=function(c) {hclust(c,method="average")}
 
mergedata.prop_rel <-merge(meta, as.matrix(data.prop_relA), by="row.names",all.x=FALSE)
 
dim(meta)

dim(data.prop_relA)


 dim(mergedata.prop_rel)


 (OTUcol1<-ncol(meta)+2)


 (OTUcol2<-ncol(mergedata.prop_rel))


 justOTU<-mergedata.prop_rel[,OTUcol1:OTUcol2]

 justOTU[1:5,1:5]

 
 rownames(justOTU[1:6,])

 rownames(justOTU)<-mergedata.prop_rel$Row.names
 
rownames(justOTU[1:6,])

 rownames(justOTU)<-mergedata.prop_rel$Row.names
 
justOTU2<-as.matrix(t(justOTU))
 
justOTU2[1:6,1:6]
      
 library(RColorBrewer)
 
HEATMAP_COl <- colorRampPalette(c('#0571b0','#92c5de','#f4a582','#ca0020'))
 
color.map2 <- function(depth) {
if (depth =="shallow") "red"
else if(depth =="deep") "darkslateblue"
 }
 
sidebarcolors <- unlist(lapply(mergedata.prop_rel$depth, color.map2))
 
clab=cbind(sidebarcolors)

colnames(clab)=c("depth")
 
pdf(file= "heatmap-sig-dif-depth.pdf",width = 25, height = 15)
 
par(mar=c(0, 2, 2, 3))
 
heatmap.3(justOTU2,  hclustfun=myclust, distfun=mydist, dendrogram="row", Rowv = TRUE, Colv = FALSE, scale="none", col = HEATMAP_COl, lhei = c(1.05, 7), margins = c(20, 40), trace="none", density.info = "none", ColSideColors= clab, keysize=1, cexRow=2.1, cexCol=2.5, cex=5.5, sepwidth=c(0.00009,0.00009), sepcolor="midnightblue", colsep=1:ncol(justOTU2), rowsep=1:nrow(justOTU2), ColSideColorsSize=4, notecex=2.5, KeyValueName="Abundance")

panel.last= legend("topright", legend=c("shallow", "deep"), fill=c("red","darkslateblue"), border=FALSE, bty="n", x.intersp = 0.4, cex=2.5, title= "depth")

dev.off()
