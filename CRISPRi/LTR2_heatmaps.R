setwd('~/Desktop/AML/AML_scripts/CRISPRi')
library(ComplexHeatmap)
library(circlize)


##heatmap plot function
##note: homer returns only copies with at least one tag, so empty copies need to be added back

draw.heatmap = function(hm.file,ltr,ltr.order=NA,sat=20,jpeg=FALSE) {
	
	##read data	
	hm = read.delim(paste('./heatmaps',hm.file,sep='/'),as.is=T)
	
	##add empty lines
	ltr.all = ltr$V4
	no.data = ltr.all[!(ltr.all %in% hm[,1])]
	hm.nam = c(hm[,1],no.data)
	hm.val = rbind(as.matrix(hm[,-1]),
	 matrix(rep(0,length(no.data)*(ncol(hm)-1)),nrow=length(no.data)))
	
	##order LTR elements
	if (is.na(ltr.order[1])) {
		total = rowSums(hm.val,)
		ltr.order = hm.nam[order(total,decreasing=T)]
	}
	mat = hm.val[match(ltr.order,hm.nam),]
	
	##define heatmap
	col = colorRamp2(c(0,sat),c('white','blue'))
	hm.draw = Heatmap(mat,cluster_rows=F,cluster_columns=F,
		 show_row_names=F,show_column_names=F,show_heatmap_legend=F,col=col)

	##draw heatmap
	if (jpeg){
		jpeg(gsub('txt','jpg',hm.file),w=1.5,h=5,units='in',res=300)
	} else {
		quartz(w=1.5,h=5)
	}
	draw(hm.draw)
	if (jpeg) dev.off()
	
	return(ltr.order)
}


##LTR annotations

ltr2b = read.delim(paste('../Annotations/LTR2B_hg38_id.bed',sep=''),header=F,as.is=T)
ltr2 = read.delim(paste('../Annotations/LTR2_hg38_id.bed',sep=''),header=F,as.is=T)


##make heatmaps

order.2b = draw.heatmap('CAS9_K562_sep_LTR2B-heatmap.txt',ltr2b)
x = draw.heatmap('CAS9_K562_LTR2B-heatmap.txt',ltr2b,ltr.order=order.2b)

order.2 = draw.heatmap('CAS9_K562_sep_LTR2-heatmap.txt',ltr2)
x = draw.heatmap('CAS9_K562_LTR2-heatmap.txt',ltr2,ltr.order=order.2)


