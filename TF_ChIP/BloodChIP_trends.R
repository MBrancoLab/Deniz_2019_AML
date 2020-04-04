##draw trend plots across LTRs
##data generated using HOMER
##uses a running mean to smooth data
##normalises signal to IgG


setwd('~/Deniz_2019_AML/BloodChIP')


##read data

flist = list.files('./bloodChIP_trends',full.names=TRUE)

data=list()
for (i in 1:length(flist)) {
	data[[i]] = read.delim(flist[i])
}


##smooth data

rmean = function(x,n) {
	cx = c(0,cumsum(x))
	rsum = (cx[(n+1):length(cx)] - cx[1:(length(cx)-n)])/n
	return(rsum)
}

smooth10 = lapply(data,function(x) rmean(x[,2],10))
pos10 = rmean(data[[i]][,1],10)


##plot profiles

par(mfrow=c(6,7),mar=rep(0.5,4))

for (ltr in c('LTR2B','LTR2C','LTR5_Hs','LTR5B','LTR12C','LTR13A')) {
	
	sub = smooth10[grep(ltr,flist)]
	mat = matrix(unlist(sub),ncol=length(sub))
	norm = mat[,-4] / mat[,4]

	for (i in 1:ncol(norm)) {
		plot(pos10,norm[,i],
		 type='l',lwd=1.5,
		 xlim=c(-1500,1500),ylim=c(0.5,3.8),
		 xlab='',ylab='',
		 xaxt='n',yaxt='n')
		abline(h=1,lty=2)
	}
}

