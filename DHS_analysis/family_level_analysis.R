##Select repeat families of interest based on enriched overlap with DHS regions
##Display enrichment of selected families across several cell types


setwd('~/Deniz_2019_AML/DHS_analysis')
library(gplots)


##read data

flist = list.files('./permutation_results/',pattern='overlap.txt')

data = list()
for (i in 1:length(flist)){
	sub = read.table(paste('./permutation_results',flist[i],sep='/'))
	sub$log.enrich[sub$log.enrich==-Inf] = -6
	sub$log.enrich[sub$log.enrich==Inf] = 6
	sub$log.enrich[is.na(sub$log.enrich)] = 0
	simple.low = grepl(')n',rownames(sub)) | grepl('rich',rownames(sub))
	data[[i]] = sub[!simple.low,]  ##exclude simple and low complexity repeats
}


##get significantly enriched repeat families

sig.enriched = lapply(data,function(x) x$pval<0.05 & x$log.enrich>1 & x$real.overlap>20)
sig.enriched = matrix(unlist(sig.enriched),ncol=length(flist))
rownames(sig.enriched) = rownames(data[[1]])
samples = unlist(lapply(strsplit(flist,'[\\._]'),function(x) x[1]))
colnames(sig.enriched) = samples


##get log2 enrichment values

log.enrich = lapply(data,function(x) x$log.enrich)
log.enrich = matrix(unlist(log.enrich),ncol=length(flist))
rownames(log.enrich) = rownames(sig.enriched)
colnames(log.enrich) = colnames(sig.enriched)


##select repeat families with enrichment in at least one cell line and >10% of AML samples

col.lines = which(colnames(sig.enriched)=='HL60' | colnames(sig.enriched)=='OCI3-1' | colnames(sig.enriched)=='Molm13-1')
sig.lines = apply(sig.enriched[,col.lines],1,any)

col.aml = grep('^S0',colnames(sig.enriched))
sig.aml = apply(sig.enriched[,col.aml],1,function(x) sum(x)>0.1*length(x))

sig.both = sig.lines & sig.aml


##display enrichments on heatmap

enr.both = log.enrich[sig.both,]

col.cd34 = grep('CD34',colnames(sig.enriched))
col.mono = grep('Monocytes',colnames(sig.enriched))
col.macro = grep('Macrophages',colnames(sig.enriched))
display.order = c(col.lines,col.aml,col.cd34,col.mono,col.macro)

heatmap.2(enr.both[,display.order],trace='none',scale='none',
 Colv=FALSE,dendrogram='row',
 col=colorRampPalette(c('blue','black','yellow'))(100),
 breaks=seq(-1,3,length.out=101),symkey=F,
 density.info='none',key.title='',key.xlab='log2 obs/exp',
 margins=c(6,6))


#############################################################


##Assi et al data

flist2 = list.files('./permutation_results/Assi_etal/',pattern='overlap.txt')

data2 = list()
for (i in 1:length(flist2)){
	sub = read.table(paste('./permutation_results/Assi_etal',flist2[i],sep='/'))
	sub$log.enrich[sub$log.enrich==-Inf] = -6
	sub$log.enrich[sub$log.enrich==Inf] = 6
	sub$log.enrich[is.na(sub$log.enrich)] = 0
	simple.low = grepl(')n',rownames(sub)) | grepl('rich',rownames(sub))
	data2[[i]] = sub[!simple.low,]  ##exclude simple and low complexity repeats
}

sig.enriched2 = lapply(data2,function(x) x$pval<0.05 & x$log.enrich>1 & x$real.overlap>20)
sig.enriched2 = matrix(unlist(sig.enriched2),ncol=length(flist2))
rownames(sig.enriched2) = rownames(data2[[1]])
samples2 = unlist(lapply(strsplit(flist2,'[\\._]'),function(x) x[1]))
colnames(sig.enriched2) = samples2

log.enrich2 = lapply(data2,function(x) x$log.enrich)
log.enrich2 = matrix(unlist(log.enrich2),ncol=length(flist2))
rownames(log.enrich2) = rownames(sig.enriched2)
colnames(log.enrich2) = colnames(sig.enriched2)

sig.aml2 = apply(sig.enriched2,1,function(x) sum(x)>0.1*length(x))
sig.both2 = sig.aml2 & names(sig.aml2) %in% names(sig.lines)[sig.lines]

enr.assi = log.enrich2[sig.both2,]

rep.order = match(c('MER57E3','U1','MSR1','LTR13','TAR1','HERVK-int','LTR2C','LTR2B','LTR13A','LTR5_Hs','LTR5B','LTR12C'),rownames(enr.assi))
 
heatmap.2(rbind(enr.assi[rep.order,],rep(NA,ncol(enr.assi)),enr.assi[-rep.order,]),
 trace='none',scale='none',
 Colv=FALSE,Rowv=FALSE,dendrogram='none',
 col=colorRampPalette(c('blue','black','yellow'))(100),
 breaks=seq(-2,4,length.out=101),symkey=F,
 density.info='none',key.title='',key.xlab='log2 obs/exp',
 margins=c(6,6))

