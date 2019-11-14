setwd('~/Desktop/AML/AML_scripts/CRISPRi')


##log2 RPM counts at dCas9 peaks

rpm = read.delim('ChIPseq_RPM.txt',as.is=T)


##Cas9 peak annotation

cas9 = read.delim('dCas9_peaks_anno.txt',as.is=T,header=F)
is.ltr2 = rpm$peak %in% cas9$V4[grep('LTR2',cas9$V11)]


##plot function

diff.plot = function(x1,x2,y1,y2) {	
	dx = x2-x1
	dy = y2-y1
	
	plot(dx,dy,pch=19,cex=0.5,las=1,
	 xlim=c(min(dx)*1.1,max(dx)*1.1),
	 ylim=c(min(dy)*1.1,max(dy)*1.1),
	 xlab='log2 FC H3K27ac',
	 ylab='log2 FC H3K9me3')

	points(dx[is.ltr2],dy[is.ltr2],pch=19,cex=0.5,col='orange')

	abline(h=0,lty=2)
	abline(v=0,lty=2)
}


##K562 plot

diff.plot(rpm$H3K27ac_K562_ctrl,
 rowMeans(cbind(rpm$H3K27ac_K562_all,rpm$H3K27ac_K562_sep)),
 rpm$H3K9me3_K562_ctrl,
 rpm$H3K9me3_K562_all)


##OCI-AML3 plot

diff.plot(rpm$H3K27ac_OCIAML3_ctrl,
 rowMeans(cbind(rpm$H3K27ac_OCIAML3_all,rpm$H3K27ac_OCIAML3_sep)),
 rpm$H3K9me3_OCIAML3_ctrl,
 rowMeans(cbind(rpm$H3K9me3_OCIAML3_all,rpm$H3K9me3_OCIAML3_sep)))







 
 