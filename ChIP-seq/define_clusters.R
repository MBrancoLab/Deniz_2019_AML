setwd('~/Desktop/AML/AML_scripts/ChIP-seq')
library(gplots)


##get overlap data

k27ac = read.delim('./LTR_histone_overlaps/allLTR_H3K27ac_overlaps.txt',as.is=T)
k4me1 = read.delim('./LTR_histone_overlaps/allLTR_H3K4me1_overlaps.txt',as.is=T)
k4me3 = read.delim('./LTR_histone_overlaps/allLTR_H3K4me3_overlaps.txt',as.is=T)
k9me3 = read.delim('./LTR_histone_overlaps/allLTR_H3K9me3_overlaps.txt',as.is=T)
overlap = list(k27ac,k4me1,k4me3,k9me3)
names(overlap) = c('H3K27ac','H3K4me1','H3K4me3','H3K9me3')


##get proportion of AML or Diff samples with mark

aml.col = grepl('AML',colnames(k27ac))
diff.col = grepl('Macro',colnames(k27ac)) | grepl('Mono',colnames(k27ac))

prop = matrix(nrow=nrow(k27ac),ncol=length(overlap)*2)
rownames(prop) = k27ac$name
colnames(prop) = rep('empty',ncol(prop))
	
for (i in 1:length(overlap)) {
	prop[,i*2-1] = rowSums(overlap[[i]][,aml.col]>0)/sum(aml.col)
	prop[,i*2] = rowSums(overlap[[i]][,diff.col]>0)/sum(diff.col)
	colnames(prop)[i*2-1] = paste(names(overlap)[i],'AML',sep='.')
	colnames(prop)[i*2] = paste(names(overlap)[i],'Diff',sep='.')
}


##heatmap plot function

plot.hist = function(data,ltr,n.clust=7,clust.order=1:n.clust) {

	##select LTR family and exclude elements with low signal
	ltr.data = data[grepl(ltr,rownames(data)) & rowSums(data)>0.05,]
	
	##k means clustering
	set.seed(10)
	k = kmeans(ltr.data,n.clust)
	
	##define order
	sorder = match(k$cluster,clust.order)
	
	##plot
	quartz(w=5,h=5)
	heatmap.2(ltr.data[order(sorder),],
	 Colv=FALSE,Rowv=FALSE,dendrogram='none',
	 trace='none',labRow=FALSE,
	 col=colorRampPalette(c('white','red'))(100),
	 breaks=seq(0,1,0.01),
	 RowSideColors=rainbow(n.clust)[sort(sorder)],
	 cexCol=0.6,key.title='')
	 
	 ##return clusters
	 return(k$cluster)
}


##plot heatmaps

ltr2b = plot.hist(prop,ltr='LTR2B',clust.order=c(4,5,1,7,3,2,6))
ltr2c = plot.hist(prop,ltr='LTR2C',n.clust=5,clust.order=c(4,5,3,1,2))
ltr5b = plot.hist(prop,ltr='LTR5B',n.clust=6,clust.order=c(1,6,5,4,3,2))
ltr5hs = plot.hist(prop,ltr='LTR5_Hs',n.clust=7,clust.order=c(3,2,5,1,7,6,4))
ltr13a = plot.hist(prop,ltr='LTR13A',n.clust=7,clust.order=c(5,1,4,7,2,3,6))
ltr12c = plot.hist(prop,ltr='LTR12C',n.clust=6,clust.order=c(6,3,5,1,4,2))


##cluster writing function

clust.write = function(ltr,clusters) {
	anno = read.delim(paste('../Annotations/',ltr,'_hg38_3kb.bed',sep=''),header=F)
	for (k in 1:max(clusters)) {
		anno.k = anno[anno$V4 %in% names(clusters)[clusters==k],]
		write.table(anno.k,
		 paste('../Annotations/',ltr,'_cluster',k,'_3kb.bed',sep=''),
		 sep='\t',quote=F,row.names=F,col.names=F)
	}
}


##write clusters

clust.write('LTR2B',ltr2b)


