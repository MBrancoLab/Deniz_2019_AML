
##Function to plot genotypes in a given order

plot.genotypes = function(sample.order) {
	
	##genotype data
	geno = read.delim('genotypes.txt',as.is=T)
	
	##reorder
	mat = as.matrix(geno[match(sample.order,gsub('-','.',geno$sampleID)),-1])

	##start plot	
	plot(NA,NA,xlim=c(0,ncol(mat)+1),ylim=c(0,nrow(mat)+1))
		
	##plot data source
	dcol = matrix(nrow=nrow(mat),ncol=2)
	dcol[mat[,1:2]] = 'lightblue'
	dcol[!mat[,1:2]] = 'white'
	
	rect(rep(0:1,each=nrow(mat)),
	 rep(nrow(mat):1,2)-1,
	 rep(1:2,each=nrow(mat)),
	 rep(nrow(mat):1,2),
	 col=dcol)
	
	
	##plot genotypes
	gcol = matrix(nrow=nrow(mat),ncol=ncol(mat)-2)
	gcol[mat[,3:ncol(mat)]] = 'orange'
	gcol[!mat[,3:ncol(mat)]] = 'white'
	gcol[is.na(mat[,3:ncol(mat)])] = 'grey'
	
	rect(rep(3:ncol(mat),each=nrow(mat)),
	 rep(nrow(mat):1,2)-1,
	 rep(3:ncol(mat),each=nrow(mat))+1,
	 rep(nrow(mat):1,2),
	 col=gcol)	
	 
}

