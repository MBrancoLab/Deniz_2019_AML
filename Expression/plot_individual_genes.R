setwd('~/Desktop/AML/AML_scripts/Expression')


##expression data

fpkm = read.delim('BP_RNA_FPKM.txt.gz',as.is=T)


##plot function

plot.gene = function(gene, log=TRUE) {
	id = which(fpkm$gene_name==gene & !is.na(fpkm$gene_name))
	sub = as.numeric(fpkm[id,-c(1,ncol(fpkm))])
	hsc = grep('stem',colnames(fpkm))-1
	mono = grep('monocyte',colnames(fpkm))-1
	macro = grep('macrophage',colnames(fpkm))-1
	aml = grep('Leukemia',colnames(fpkm))-1
	
	par(mar=c(6,4,2,2))
	if (log) {
		sub = log2(sub+0.01)
		las = 2
		ylab = 'Log2 Expression'
	} else {
		las = 3
		ylab = 'Expression'
	}
	boxplot(sub[hsc],sub[mono],sub[macro],sub[aml],las=las,col='wheat',
	 names=c('HSC','Monocyte','Macrophage','AML'),ylab=ylab,
	 lty=1,pch=19,cex=0.5,main=gene)
}


##plot genes of interest

plot.gene('RPS14')



