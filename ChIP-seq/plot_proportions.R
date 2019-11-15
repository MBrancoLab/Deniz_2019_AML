##Plots percentage of elements from each LTR family that overlap different histone marks

setwd('~/Deniz_2019_AML/ChIP-seq')


##get overlap data

k27ac = read.delim('./LTR_histone_overlaps/allLTR_H3K27ac_overlaps.txt',as.is=T)
k4me1 = read.delim('./LTR_histone_overlaps/allLTR_H3K4me1_overlaps.txt',as.is=T)
k4me3 = read.delim('./LTR_histone_overlaps/allLTR_H3K4me3_overlaps.txt',as.is=T)
k9me3 = read.delim('./LTR_histone_overlaps/allLTR_H3K9me3_overlaps.txt',as.is=T)
overlap = list(k27ac,k4me1,k4me3,k9me3)
names(overlap) = c('H3K27ac','H3K4me1','H3K4me3','H3K9me3')


##simplify

overmat = lapply(overlap, function(x) {
	mat = as.matrix(x[,-c(1:6)]>0)
	rownames(mat) = unlist(lapply(strsplit(x$name,split='_'),function(y) y[1]))
	return(mat)
})


##add promoter/enhancer classification

overmat$promoter = overmat$H3K27ac & overmat$H3K4me3
overmat$enhancer = overmat$H3K27ac & overmat$H3K4me1


##get % elements per family

prop = lapply(overmat, function(x) {
	fam = factor(rownames(x))
	p = apply(x,2,function(y) {
		tapply(y,fam,sum)/tapply(y,fam,length)*100
	})
})



##plot % of elements with mark

mark = 'H3K4me1'

mdata = prop[[mark]]
aml.col = grepl('AML',colnames(mdata))
diff.col = grepl('Macro',colnames(mdata)) | grepl('Mono',colnames(mdata))

quartz(w=5.2,h=3.8)
plot(NA,NA,xlim=c(0.5,12.5),ylim=c(0,42),las=1,xaxt='n',
 xlab='',ylab='% of elements',main=mark)
for (i in 1:nrow(mdata)) {
	points(rnorm(sum(aml.col),mean=i*2-1,sd=0.1),mdata[i,aml.col],pch=19,cex=0.5,col='black')
	lines(c(i*2-1.2,i*2-0.8),rep(mean(mdata[i,aml.col]),2),lwd=2,col='blue')
	points(rnorm(sum(diff.col),mean=i*2,sd=0.1),mdata[i,diff.col],pch=19,cex=0.5,col='grey')
	lines(c(i*2-0.2,i*2+0.2),rep(mean(mdata[i,diff.col]),2),lwd=2,col='blue')
}
axis(1,at=1:12,labels=rep(c('AML','Diff.'),6),las=2)


##Wilcoxon tests

p = apply(mdata,1,function(x) wilcox.test(x[aml.col],x[diff.col])$p.val)
p.adjust(p,method='BH')


