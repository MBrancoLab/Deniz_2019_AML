##Makes DHS profiles across LTRs using HOMER data


setwd('~/Deniz_2019_AML/DHS_analysis')
library(ComplexHeatmap)
library(circlize)


##heatmap plot function
##note: homer returns only copies with at least one tag, so empty copies need to be added back

draw.heatmap = function(files,te,jpeg=FALSE) {
	ltr.all = read.delim(paste('../Annotations/',te,'_hg38_3kb.bed',sep=''),header=F,as.is=T)$V4
	n = length(files)
	
	hm.list = NULL
	for (i in 1:n) {
	##add empty lines
		hm = read.delim(files[i],as.is=T)
		no.data = ltr.all[!(ltr.all %in% hm[,1])]
		hm.nam = c(hm[,1],no.data)
		hm.val = rbind(as.matrix(hm[,-1]),
		 matrix(rep(0,length(no.data)*(ncol(hm)-1)),nrow=length(no.data)))
		
	##make heatmap list object
		total = rowSums(hm.val,)
		if (i==1) ltr.order = hm.nam[order(total,decreasing=T)]
		mat = hm.val[match(ltr.order,hm.nam),]
		col = colorRamp2(c(min(mat,na.rm=T),quantile(mat,0.97,na.rm=T)),c('white','blue'))
		hm.list = hm.list + Heatmap(mat,cluster_rows=F,cluster_columns=F,
		 show_row_names=F,show_column_names=F,show_heatmap_legend=F,col=col)
	}
	
	##plot
	if (jpeg){
		jpeg(paste(te,'heatmap.jpg',sep='_'),w=1.6,h=5,units='in',res=300)
	} else {
		quartz(w=1.6,h=5)
	}
	draw(hm.list)
	if (jpeg) dev.off()
	return(ltr.order)
}


##HOMER output files (provided as zip file in the repository)

flist = list.files('./heatmaps',full.names=TRUE)


##generate heatmap profiles

draw.heatmap(flist[grep('LTR2B-OCI3',flist)],'LTR2B',jpeg=T)
draw.heatmap(flist[grep('LTR2C-OCI3',flist)],'LTR2C',jpeg=T)
draw.heatmap(flist[grep('LTR5B-OCI3',flist)],'LTR5B',jpeg=T)
draw.heatmap(flist[grep('LTR5_Hs-OCI3',flist)],'LTR5_Hs',jpeg=T)
draw.heatmap(flist[grep('LTR12C-OCI3',flist)],'LTR12C',jpeg=T)
draw.heatmap(flist[grep('LTR13A-OCI3',flist)],'LTR13A',jpeg=T)



