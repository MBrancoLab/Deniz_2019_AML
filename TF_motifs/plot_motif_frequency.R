##Takes output from FIMO to calculate motif frequency in each LTR family
##Plots frequency of selected motifs and how it compares between DHS+ and DHS- copies


setwd('~/Deniz_2019_AML/TF_motifs')


##read in fimo files

fimo.files = list.files('./fimo_files')

fimo = list()
for (i in 1:length(fimo.files)) {
	fimo[[i]] = read.delim(paste('./fimo_files',fimo.files[i],sep='/'),as.is=T)
}
names(fimo) = gsub('_fimo.txt','',fimo.files)


##enriched motifs table

enr = read.table('enriched_motifs.txt',row.names=1)


##DHS elements

dhs = read.delim('../DHS_analysis/allLTR_DHS_overlaps.txt',as.is=T)
is.dhs = rowSums(dhs[,7:ncol(dhs)])>=5


##get motif frequency

freq = lapply(names(fimo), function(x) {
	f.non = numeric(nrow(enr))
	f.dhs = numeric(nrow(enr))
	for (i in 1:nrow(enr)) {
		motif.elements = fimo[[x]]$sequence.name[fimo[[x]]$X.pattern.name==rownames(enr)[i]]
		f.non[i] = sum(unique(motif.elements) %in% dhs$name[!is.dhs])/sum(grepl(x,dhs$name[!is.dhs]))
		f.dhs[i] = sum(unique(motif.elements) %in% dhs$name[is.dhs])/sum(grepl(x,dhs$name[is.dhs]))
	}
	res = cbind(f.non,f.dhs)
	rownames(res) = rownames(enr)
	colnames(res) = c('nonDHS','DHS')
	return(res)
})
names(freq) = names(fimo)


##scatter plots of dhs vs non-dhs

dhs.plot = function(ltr) {
	x = freq[[ltr]]*100
	plot(x[,1],x[,2],pch=19,cex=0.5,las=1,
	 xlim=c(0,100),ylim=c(0,100),
	 xlab='% non DHS elements',ylab='% DHS elements',
	 main=ltr,col='grey')
	abline(0,1,lty=2)
	
	##option to highlight enriched motifs from AME analysis:
	#sig = enr[,colnames(enr)==paste(ltr,'dhs',sep='_')]
	#points(x[sig,1],x[sig,2],pch=19,cex=0.5,col='red')
	#text(x[sig,1],x[sig,2],rownames(enr)[sig],cex=0.5,pos=4)
	
	##option to highlight specific motifs:
	#sig=rownames(enr) %in% c('CEBPB_HUMAN.H11MO.0.A','GATA2_HUMAN.H11MO.0.A','HXA9_HUMAN.H11MO.0.B',
	# 'MEIS1_HUMAN.H11MO.0.A','SPI1_HUMAN.H11MO.0.A','TAL1_HUMAN.H11MO.0.A')
	#points(x[sig,1],x[sig,2],pch=19,cex=0.7,col='red')
	#return(x[sig,])
}

quartz(w=4.5,h=5)
dhs.plot('LTR2B')
dhs.plot('LTR2C')
dhs.plot('LTR5B')
dhs.plot('LTR5_Hs')
dhs.plot('LTR12C')
dhs.plot('LTR13A')


##barplots of motifs of interest

moi = c('SPI1_HUMAN.H11MO.0.A','CEBPB_HUMAN.H11MO.0.A','GATA2_HUMAN.H11MO.0.A','HXA9_HUMAN.H11MO.0.B',
 'MEIS1_HUMAN.H11MO.0.A','IKZF1_HUMAN.H11MO.0.C','NFYA_HUMAN.H11MO.0.A','PKNX1_HUMAN.H11MO.0.B',
 'RUNX1_HUMAN.H11MO.0.A','STA5A_HUMAN.H11MO.0.A','TAL1_HUMAN.H11MO.0.A')

freq.sel = lapply(names(fimo), function(x) {
	f = numeric(length(moi))
	for (i in 1:length(moi)) {
		motif.elements = fimo[[x]]$sequence.name[fimo[[x]]$X.pattern.name==moi[i]]
		f[i] = length(unique(motif.elements))/sum(grepl(x,dhs$name))
	}
	return(f)
})

fmat = matrix(unlist(freq.sel),ncol=length(freq.sel))
rownames(fmat) = gsub('_HUMAN.H11MO[[:print:]]+','',moi)
colnames(fmat) = names(fimo)

barplot(t(fmat)*100,beside=T,las=2,ylab='% of elements',
 col=c('orange','wheat','red3','tomato','navy','lightblue'))

