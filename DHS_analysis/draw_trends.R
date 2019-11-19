##draw trend plots across LTRs
##data generated using HOMER
##uses a running mean to smooth data


setwd('~/Deniz_2019_AML/DHS_analysis')


##read data

flist = list.files('./trends',full.names=TRUE)

data=list()
for (i in 1:length(flist)) {
	data[[i]] = read.delim(flist[i])
}


##define sample groups

clines = grepl('HL60',flist) | grepl('Molm13',flist) | grepl('OCI3',flist) 
macro = grepl('C005VG45',flist) | grepl('S001S745',flist) | grepl('S0022I44',flist)
mono = grepl('C0010K46',flist) | grepl('C001UY46',flist)
aml = !clines & !macro & !mono


##running mean function

rmean = function(x,n) {
	cx = c(0,cumsum(x))
	rsum = (cx[(n+1):length(cx)] - cx[1:(length(cx)-n)])/n
	return(rsum)
}



##smooth data

smooth50 = lapply(data,function(x) rmean(x[,2],50))
smooth10 = lapply(data,function(x) rmean(x[,2],10))


##LTR to plot

ltr = grepl('LTR2B',flist)


##plot Blueprint profiles

aml.cov = matrix(unlist(smooth50[aml&ltr]),ncol=sum(aml&ltr))
macro.cov = matrix(unlist(smooth50[macro&ltr]),ncol=sum(macro&ltr))
mono.cov = matrix(unlist(smooth50[mono&ltr]),ncol=sum(mono&ltr))

pos = rmean(data[[i]][,1],50)

quartz(w=3.5,h=4)
plot(NA,NA,xlim=c(-1000,1000),ylim=c(0.3,1.8),
 xlab='Position from LTR centre (bp)',
 ylab='Relative DNAse-seq coverage')
for (i in 1:ncol(aml.cov)) lines(pos,aml.cov[,i],col='grey')
for (i in 1:ncol(macro.cov)) lines(pos,macro.cov[,i],col='blue')
for (i in 1:ncol(mono.cov)) lines(pos,mono.cov[,i],col='red')


##plot cell line profiles

clines.cov = matrix(unlist(smooth10[clines&ltr]),ncol=sum(clines&ltr))

pos2 = rmean(data[[i]][,1],10)

quartz(w=3.5,h=4)
plot(NA,NA,xlim=c(-1000,1000),ylim=c(0.1,0.7),
 xlab='',
 ylab='Relative DNAse-seq coverage')
cols = c('black','orange','lightblue')
for (i in 1:ncol(clines.cov)) lines(pos2,clines.cov[,i],col=cols[i],lwd=1.5)


