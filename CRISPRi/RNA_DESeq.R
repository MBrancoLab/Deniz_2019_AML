setwd('~/Desktop/AML/AML_scripts/CRISPRi')
library('DESeq2')

##Data

reads = read.table('LTR2B_RNA_raw.txt', header = T, row.names = 1)
colData = DataFrame(sample=factor(rep(c('K562.ctrl','K562.krab','OCI.ctrl','OCI.krab'),c(2,3,2,3))))
colData1 = DataFrame(sample=factor(rep(c('Ctrl','KRAB'),c(2,3))))


##Expression values (VST)

##dds = DESeqDataSetFromMatrix(reads,colData,design=~sample)
##vsd = varianceStabilizingTransformation(dds)
##expr = assay(vsd)
expr = read.table('LTR2B_RNA_vsd.txt',as.is=T,row.names=1)


##Get DE genes

dds1 = DESeqDataSetFromMatrix(reads[,1:5],colData1,design=~sample)
dds1 = DESeq(dds1)
res1 = results(dds1,contrast=c('sample','KRAB','Ctrl'))
k562 = as.data.frame(subset(res1,padj<0.05))
k562 = k562[order(k562$log2FoldChange,decreasing=T),]

dds2 = DESeqDataSetFromMatrix(reads[,6:10],colData1,design=~sample)
dds2 = DESeq(dds2)
res2 = results(dds2,contrast=c('sample','KRAB','Ctrl'))
oci = as.data.frame(subset(res2,padj<0.05))
oci = oci[order(oci$log2FoldChange,decreasing=T),]

k562 = read.delim('LTR2B_de_k562.txt', row.names=1)
oci = read.delim('LTR2B_de_oci.txt', row.names=1)



##Check expression

plot.gene = function(gene,log=TRUE) {
	id = rownames(expr)==gene
	sample = colData$sample
	
	if (log) {
		data = expr[id,]
	} else {
		data = 2^expr[id,]
	}
	
	av = tapply(as.numeric(data),sample,mean)
	h=barplot(av,main=gene,ylab='Relative expression',col='grey',ylim=c(0,max(data)*1.1),las=2)
	for (i in 1:nlevels(sample)) {
		samp = sample==levels(sample)[i]
		points(rep(h[i],sum(samp)),data[samp],pch=19,cex=0.5)
	}
	return(data)
}

plot.gene('ZNF320')
plot.gene('APOC1')
plot.gene('APOC2')
plot.gene('APOC4-APOC2')
plot.gene('APOE')
plot.gene('APOL1')
plot.gene('IL23R')
