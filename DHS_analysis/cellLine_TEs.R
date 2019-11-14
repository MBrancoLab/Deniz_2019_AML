

##get cell line data

data = list(hl60=read.table('./permutation_results/HL60_filt.npf_overlap.txt'),
 molm13=read.table('./permutation_results/Molm13-1_filt.npf_overlap.txt'),
 ociaml3=read.table('./permutation_results/OCI3-1_filt.npf_overlap.txt'))


##get numbers for families of interest

ltr = c('LTR2B','LTR2C','LTR5B','LTR5_Hs','LTR13A','LTR12C')

subdata = lapply(data,function(x) x[match(ltr,rownames(x)),])
abs.n = lapply(subdata,function(x) x$real.overlap)
rel.n = lapply(subdata,function(x) x$real.overlap/x$total*100)

abs.mat = matrix(unlist(abs.n),nrow=length(abs.n),byrow=T)
rel.mat = matrix(unlist(rel.n),nrow=length(rel.n),byrow=T)


##plots

quartz(w=4.5,h=4.5)
h = barplot(rel.mat,beside=T,names.arg=ltr,las=2,ylab='% of elements in family',
 ylim=c(0,32),col=c('black','orange','lightblue'))
text(h,rel.mat+1.7,labels=abs.mat,srt=90,cex=0.8)
legend(2,30,legend=c('HL60','Molm13','OCI-AML3'),fill=c('black','orange','lightblue'),bty='n')

h = barplot(abs.mat,beside=T,names.arg=ltr,las=2,ylab='n elements',
 col=c('black','orange','lightblue'))
legend(2,350,legend=c('HL60','Molm13','OCI-AML3'),fill=c('black','orange','lightblue'),bty='n')

