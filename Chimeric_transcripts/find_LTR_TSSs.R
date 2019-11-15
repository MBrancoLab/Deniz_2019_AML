##Finds overlaps between LTRs of interest and Stringtie TSSs
##Uses bedtools

setwd('~/Deniz_2019_AML/Chimeric_transcripts')
library('gplots')
path_to_intersectbed = '~/Documents/bedtools-2.26.0/bin/intersectBed'


##load TSSs

load('Blueprint_TSSs.Rdata')


##overlap LTRs with TSSs

overlap.tss = function(ltr.file,tss) {
	tss.bed = cbind(tss,tss$pos+1)
	b.file = tempfile()
	out.file = tempfile()
	write.table(tss.bed,b.file,sep='\t',quote=F,col.names=F,row.names=F)
	command = paste(path_to_intersectbed,'-c -a',ltr.file,'-b',b.file,'>',out.file)
	try(system(command))
	out = read.delim(out.file,header=F,as.is=T)
	unlink(c(b.file,out.file))
	return(out$V7>0)
}

ltr.file = '../Annotations/all_LTRs.bed'
ltr = lapply(tss,function(x) overlap.tss(ltr.file,x))


##count putative chimeric transcripts per family

ltr.bed = read.delim(ltr.file,header=F,as.is=T)

ltr2b.sum = unlist(lapply(ltr,function(x) sum(x[grep('LTR2B',ltr.bed$V4)])))
ltr2c.sum = unlist(lapply(ltr,function(x) sum(x[grep('LTR2C',ltr.bed$V4)])))
ltr5hs.sum = unlist(lapply(ltr,function(x) sum(x[grep('LTR5_Hs',ltr.bed$V4)])))
ltr5b.sum = unlist(lapply(ltr,function(x) sum(x[grep('LTR5B',ltr.bed$V4)])))
ltr12c.sum = unlist(lapply(ltr,function(x) sum(x[grep('LTR12C',ltr.bed$V4)])))
ltr13a.sum = unlist(lapply(ltr,function(x) sum(x[grep('LTR13A',ltr.bed$V4)])))

tss.sum = rbind(ltr2b.sum,ltr2c.sum,ltr5hs.sum,ltr5b.sum,ltr12c.sum,ltr13a.sum)
colnames(tss.sum) = unlist(lapply(strsplit(names(ltr),'_'),function(x) x[1]))

sel.cols = colnames(tss.sum)=='AML' | colnames(tss.sum)=='Macro' | colnames(tss.sum)=='Mono'

barplot(tss.sum[,sel.cols],las=2,cex.names=0.8,ylab='n TSSs at LTRs',
 col=c('black','grey','orange','tomato','lightblue','wheat'))


##select AML-associated transcripts 

ltr.mat = matrix(as.numeric(unlist(ltr)),ncol=length(ltr))
colnames(ltr.mat) = colnames(tss.sum)
rownames(ltr.mat) = ltr.bed$V4

in.diff = rowSums(ltr.mat[,colnames(ltr.mat)=='Macro'|colnames(ltr.mat)=='Mono'])>=1
in.aml = rowSums(ltr.mat[,colnames(ltr.mat)=='AML'])>=2
aml.only = in.aml & !in.diff

aml.ltr = cbind(ltr.bed[aml.only,],aml.count=rowSums(ltr.mat[aml.only,colnames(ltr.mat)=='AML']))
write.table(aml.ltr,'AML_associated_TSS_LTRs.txt',sep='\t',quote=F,row.names=F)

