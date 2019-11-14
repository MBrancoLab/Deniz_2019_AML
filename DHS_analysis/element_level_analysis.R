
setwd('~/Desktop/AML/AML_scripts/DHS_analysis')
library(gplots)
source('plot_genotypes.R')


##get overlap data

overlap = read.delim('allLTR_DHS_overlaps.txt',as.is=T)


##DHS matrix

dhs.mat = as.matrix(overlap[,-c(1:6)])
dhs.mat[dhs.mat>1] = 1
rownames(dhs.mat) = overlap$name
colnames(dhs.mat) = colnames(overlap[-c(1:6)])

line.cols = colnames(dhs.mat)=='HL60' | colnames(dhs.mat)=='Molm13.1' | colnames(dhs.mat)=='OCI3.1'
cd34.cols = grepl('CD34',colnames(dhs.mat))
aml.cols = grepl('^S0',colnames(dhs.mat))
diff.cols = grepl('Macrophages',colnames(dhs.mat)) | grepl('Monocytes',colnames(dhs.mat))
assi.cols = !line.cols & !cd34.cols & !aml.cols & !diff.cols


##plot DHS matrix

one.dhs = rowSums(dhs.mat[,aml.cols|assi.cols|diff.cols])>=1

hm = heatmap.2(dhs.mat[one.dhs,],trace='none',key=FALSE,
 col=colorRampPalette(c('black','yellow'))(2),
 dendrogram='row',labRow=NA,
 hclustfun = function(x) hclust(x,method='average'))


##plot number of DHS elements

plot(NA,NA,xlim=c(0,7),ylim=c(0,1200),las=1,xaxt='n',xlab='',ylab='n of DHS elements')
points(rnorm(sum(aml.cols),mean=1,sd=0.2),colSums(dhs.mat[,aml.cols]),pch=19,cex=0.5)
points(rnorm(sum(diff.cols),mean=3,sd=0.2),colSums(dhs.mat[,diff.cols]),pch=19,cex=0.5)
points(rnorm(sum(assi.cols),mean=5,sd=0.2),colSums(dhs.mat[,assi.cols]),pch=19,cex=0.5)
axis(1,at=c(1,3,5),labels=c('AML','Diff.','Assi'))


##correlation matrix (with genotypes)

cor.mat = cor(dhs.mat[one.dhs,aml.cols|assi.cols|diff.cols])

hm2 = heatmap.2(cor.mat,trace='none',
 col=colorRampPalette(c('blue','white','red'))(100),
 hclustfun = function(x) hclust(x,method='complete'),
 breaks=seq(0,1,0.01))

n.dhs = colSums(dhs.mat[,aml.cols|assi.cols|diff.cols])
barplot(n.dhs[rev(hm2$rowInd)],las=2,cex.names=0.6,space=0,cex.axis=0.7)

plot.genotypes(rownames(cor.mat)[rev(hm2$rowInd)])


##function for extracting correlations between samples with the same mutation

geno = read.delim('genotypes.txt',as.is=T)

geno.cor = function(cor.mat,geno,mutation) {
	
	##reorder
	gen.mat = as.matrix(geno[match(rownames(cor.mat),gsub('-','.',geno$sampleID)),-1])

	##get mutation data
	mut = which(gen.mat[,colnames(gen.mat)==mutation])

	##get correlation coefficients
	gen.cor = numeric()
	for (i in 1:(length(mut)-1)) {
		for (j in (i+1):length(mut)) {
			gen.cor = c(gen.cor,cor.mat[mut[i],mut[j]])
		}
	}
	
	##as above, for controls
	ctrl = (1:nrow(gen.mat))[-mut]
	ctrl.cor = numeric()
	for (i in 1:(length(ctrl)-1)) {
		for (j in (i+1):length(ctrl)) {
			ctrl.cor = c(ctrl.cor,cor.mat[ctrl[i], ctrl[j]])
		}
	}

	return(list(gen.cor,ctrl.cor))
}


##analyse correlations for different mutations

all.gcor = c(geno.cor(cor.mat,geno,'FLT3.ITD'),
 geno.cor(cor.mat,geno,'NPM1'),
 geno.cor(cor.mat,geno,'DNMT3A'),
 geno.cor(cor.mat,geno,'CEBPA'),
 geno.cor(cor.mat,geno,'RUNX1'),
 geno.cor(cor.mat,geno,'t8.21'))

boxplot(all.gcor,ylim=c(-0.1,0.8),
 at=1:length(all.gcor) + rep(seq(from=0,by=0.5,length.out=length(all.gcor)/2),each=2),
 lty=1,pch=19,cex=0.5,col=c('orange','grey'),
 xaxt='n',las=1,ylab='correlation coefficient')
axis(1,at=seq(from=1.5,by=2.5,length.out=length(all.gcor)/2),
 labels=c('FLT3-ITD','NPM1','DNMT3A','CEBPA','RUNX1','t(8;21)'),las=2)

p = numeric(length(all.gcor)/2)
for (i in 1:length(p)) {
	p[i] = t.test(all.gcor[[i*2-1]],all.gcor[[i*2]])$p.value
}
padj = p.adjust(p,method='BH')
