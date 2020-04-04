##Generates trend plots across LTR2B elements from each of the clusters defined in define_clusters.R

setwd('~/Deniz_2019_AML/ChIP-seq')


##data files (provided as a zip file in this repository)

flist = list.files('./trends',full.names=TRUE)


##plot function


plot.trend = function(cluster,hist) {

	##select cluster and histone mark
	fsel = flist[grepl(cluster,flist) & grepl(hist,flist)]
	
	
	##read data
	data=list()
	for (i in 1:length(fsel)) {
		data[[i]] = read.delim(fsel[i])
	}
	
	##make coverage matrix
	cov = matrix(unlist(lapply(data,function(x) x[,2])),ncol=length(data))
	cov.names = unlist(lapply(data,function(x) colnames(x)[2]))
	colnames(cov) = gsub('\\.ERX[[:print:]]+','',cov.names)
	
	
	##define cell types
	diff.col = colnames(cov) %in% c('C000S5H2','C0011IH1','C001UYH2','C00264H1','C00280H2',
	 'C004SQH1','C005PSH2','C005VGH1','S000RDH2','S001S7H1','S0022IH1','S00390H1','S00BHQH1',
	 'S00DVRH1','S01F8KH1')
	aml.col = !diff.col
	
	
	##get coverage means
	aml.av = rowMeans(cov[,aml.col])
	diff.av = rowMeans(cov[,diff.col])
	
	
	##plot
	y.lim = c(min(c(aml.av,diff.av)),max(c(aml.av,diff.av,2)))
	x = data[[1]][,1]
	plot(x,aml.av,type='l',lwd=2,ylim=y.lim,xlab='',ylab='coverage')
	lines(x,diff.av,lwd=2,col='grey')
}



##plot trends

for (i in 1:7) {
	quartz(w=3,h=3.5)
	plot.trend(paste('cluster',i,sep=''),'H3K9me3')
}


