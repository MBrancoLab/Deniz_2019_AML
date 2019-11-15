##Merges expression and DHS data
##Generates plots to check the association between DHS+ LTRs and gene expression

setwd('~/Deniz_2019_AML/Expression')
library(vioplot)


##expression data

fpkm = read.delim('BP_RNA_FPKM.txt.gz',as.is=T)


##DHS data

dhs = read.delim('~/Desktop/AML/AML_scripts/DHS_analysis/allLTR_DHS_overlaps.txt',as.is=T)



######################
#### prepare data ####



##get average expression values for cell groups

hsc = rowMeans(log2(fpkm[,grep('stem.cell',colnames(fpkm))]+0.01))
mono = rowMeans(log2(fpkm[,grep('monocyte',colnames(fpkm))]+0.01))
macro = rowMeans(log2(fpkm[,grep('macrophage',colnames(fpkm))]+0.01))
aml = rowMeans(log2(fpkm[,grep('Leukemia',colnames(fpkm))]+0.01))


##make matching dhs and expression matrices for AML samples

meta = read.delim('../blueprint_files.tsv',as.is=T)
dhs.meta = meta[meta$Experiment=='DNase-Seq' & meta$Format=='BED' & meta$Sub.group=='Acute Myeloid Leukemia',]

for (i in 1:ncol(dhs)) {
	id = grep(colnames(dhs)[i],dhs.meta$URL)[1]
	if (!is.na(id)) colnames(dhs)[i] = dhs.meta$Donor[id]
}

expr.donor = gsub('Acute.Myeloid.Leukemia.','',colnames(fpkm))

aml.dhs = as.matrix(dhs[,colnames(dhs) %in% expr.donor]>=1)
aml.expr = log2(as.matrix(fpkm[,match(colnames(aml.dhs),expr.donor)])+0.01)


##get average expression values from dhs and non-dhs AML groups

av.dhs = numeric(nrow(aml.expr))
av.nodhs = numeric(nrow(aml.expr))
for (i in 1:nrow(aml.expr)) {
	if (is.na(fpkm$LTR_name[i])) {
		av.dhs[i] = NA
		av.nodhs[i] = NA
	} else {
		is.dhs = aml.dhs[fpkm$LTR_name[i]==dhs$name]
		av.dhs[i] = mean(aml.expr[i,is.dhs])
		av.nodhs[i] = mean(aml.expr[i,!is.dhs])
	}
}


##make LTRs groups based on DHS data

diff.dhs = dhs[,grepl('Monocytes',colnames(dhs)) | grepl('Macrophages',colnames(dhs))]>=1

nodhs.ltrs = dhs$name[rowSums(aml.dhs)==0 & rowSums(diff.dhs)==0]
aml.ltrs = dhs$name[rowSums(aml.dhs)>=2 & rowSums(diff.dhs)==0]
ubi.ltrs = dhs$name[rowSums(aml.dhs)>=2 & rowSums(diff.dhs)>=2]


##set distance threshold

near = fpkm$LTR_distance <= 50000



##############################################################
#### plot average AML expression for different LTR groups ####


vioplot(aml[near & fpkm$LTR_name %in% nodhs.ltrs],
 aml[near & fpkm$LTR_name %in% aml.ltrs],
 aml[near & fpkm$LTR_name %in% ubi.ltrs],
 names=c('None','AML','Ubi'),
 col='wheat',yaxt='n',ylab='log2 FPKM')
axis(2,seq(-5,15,5),las=1)



##############################################################
#### compare expression between DHS+ and DHS- AML samples ####

plus = av.dhs[near & fpkm$LTR_name %in% aml.ltrs]
minus = av.nodhs[near & fpkm$LTR_name %in% aml.ltrs]

plot(minus,plus,pch=19,cex=0.5,col='grey',las=1,
 xlab='DHS- AMLs',ylab='DHS+ AMLs')
abline(0,1,lty=2,col='blue')


##pinpoint DHS-associated interest

sel = (plus>0 | minus>0) & plus-minus>2
points(minus[sel],plus[sel],pch=19,cex=0.5,col='red')

aml.genes = fpkm$gene_name[near & fpkm$LTR_name %in% aml.ltrs]
goi = cbind(aml.genes,minus,plus)[sel,]



#########################################################
#### plot expression of AML LTRs in other cell types ####
##(not included in the paper)

exp.list = list(hsc,mono,macro,av.dhs,av.nodhs)
aml.list = lapply(exp.list,function(x) x[near & fpkm$LTR_name %in% aml.ltrs])

vioplot(aml.list,names=c('HSC','Mono','Macro','DHS+','DHS-'),
 col='wheat',yaxt='n',ylab='log2 FPKM')
axis(2,seq(-5,15,5),las=1)


