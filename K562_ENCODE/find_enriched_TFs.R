setwd('~/Desktop/AML/AML_scripts/K562_ENCODE')
library(gplots)


##read overlap data

enc = read.table('LTR_hg19_overlaps.txt',row.names=1)
ran = read.table('LTR_shuffled_overlaps.txt',row.names=1)


##exclude datasets with treatments

meta = read.delim('K562_metadata.tsv',as.is=T)
bed = substr(rownames(enc),gregexpr('ENC',rownames(enc)),nchar(rownames(enc))-4)
meta.bed = meta[match(bed,meta$File.accession),]
untreated = meta.bed$Biosample.treatments==''

enc = enc[untreated,]
ran = ran[untreated,]


##Fisher's test function

fish.test = function(x) {
	obs = x[1]
	exp = x[2]
	total = x[3]
	mat = matrix(c(obs,exp,total-obs,total-exp),ncol=2)
	p = fisher.test(mat)$p.value
	return(p)
}


##find files with significant enrichment

n.ltr = matrix(rep(c(2740,188,326,295,431,645),nrow(enc)),byrow=T,nrow=nrow(enc))  #see 'LTR_counts.txt' file

pval = matrix(nrow=nrow(enc),ncol=ncol(enc))
for (i in 1:ncol(enc)) {
	test.mat = cbind(enc[,i],ran[,i],n.ltr[,i])
	p = apply(test.mat,1,fish.test)
	pval[,i] = p.adjust(p,method='BH')
}

sig = pval<0.05 & enc>ran & enc/n.ltr>0.05


##plot enriched files

par(mfrow=c(2,3))
for (i in 1:ncol(enc)) {
	plot(ran[,i],enc[,i],pch=19,cex=0.3,col='grey',main=colnames(enc)[i],
	 xlab='Expected',ylab='Observed',las=1)
	points(ran[sig[,i],i],enc[sig[,i],i],pch=19,cex=0.3,col='red')
	abline(0,1,lty=2)
}
#(Note: the two groups of points are from different processing pipelines)


##get average enrichment for significant TFs

is.sig = apply(sig,1,any) & !grepl('^H[[:digit:]][[:print:]]+',rownames(sig))  #exlcude histone marks
oe.sig = log2(enc[is.sig,]+1) - log2(ran[is.sig,]+1)

no.gfp = gsub('eGFP-','',rownames(oe.sig)) #remove eGFP prefix from GFP fusion experiments
tf = factor(gsub('-ENC[[:print:]]+','',no.gfp)) #217 TFs

oe.av = matrix(nrow=nlevels(tf),ncol=ncol(oe.sig))
for (i in 1:ncol(oe.sig)) oe.av[,i] = tapply(oe.sig[,i],tf,mean)
rownames(oe.av) = levels(tf)
colnames(oe.av) = colnames(oe.sig)


##add info on TF expression in AML

fpkm = read.delim('../Expression/BP_RNA_FPKM.txt.gz',as.is=T)

fpkm$gene_name[fpkm$gene_name=='NSD2'] = 'WHSC1' #special cases where gene name matching failed
fpkm$gene_name[fpkm$gene_name=='EMSY'] = 'C11orf30'
fpkm$gene_name[fpkm$gene_name=='CAVIN1'] = 'PTRF'
fpkm$gene_name[fpkm$gene_name=='AC016586.1'] = 'SIRT6'
fpkm$gene_name[fpkm$gene_name=='AC118549.1'] = 'ZZZ3'

tf.exp = fpkm[match(rownames(oe.av),fpkm$gene_name),]
aml.max = apply(tf.exp[,grep('Leukemia',colnames(fpkm))],1,max)
oe.expr = cbind(oe.av,aml.max)

write.table(oe.expr,'enriched_TFs.txt',sep='\t',quote=F,col.names=NA)


##plot selection of TFs

sel.tf = c('NFYA','NFYB','TAL1','RUNX1','CEBPB','CEBPD','GATA2','SPI1','MAX','MYC',
 'NFE2','IRF1','E2F7','RELA','ZEB2','JUN','NCOR1','IKZF1','PKNOX1','ARNT',
 'LEF1','ETV6','EGR1','BCOR','ELF1','ELF4','STAT5A','ETS2')

quartz(w=7,h=4)
heatmap.2(t(oe.av[match(sel.tf,rownames(oe.av)),]),
 trace='none',
 Rowv=NA,dendrogram='none',
 margins=c(4.5,6.5),cexCol=0.8,
 col=colorRampPalette(c('blue','white','red'))(100),
 breaks=seq(-4,4,length.out=101))


