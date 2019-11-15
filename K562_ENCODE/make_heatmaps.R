##Makes TF ChIP-seq profiles from HOMER data


setwd('~/Deniz_2019_AML/K562_ENCODE')
library(ComplexHeatmap)
library(circlize)


##heatmap plot function
##note: homer returns only copies with at least one tag, so empty copies need to be added back

draw.heatmap = function(hm.file,ltr,ltr.order=NA,jpeg=FALSE) {
	
	##read data	
	hm = read.delim(paste('./LTR_heatmaps',hm.file,sep='/'),as.is=T)
	
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
	col = colorRamp2(c(0,4),c('white','blue'))   #define saturation point here
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


##make heatmaps

fam.hm = function(fam,order.by) {
	ltr = read.delim(paste('../Annotations/',fam,'_hg38_id.bed',sep=''),header=F,as.is=T)
	hm.files = list.files('./LTR_heatmaps',pattern=fam)
	
	order.file = grep(order.by,hm.files)
	ltr.order = draw.heatmap(hm.files[order.file],ltr,jpeg=TRUE)
	
	for (f in hm.files[-order.file]) {
			temp = draw.heatmap(f,ltr,ltr.order=ltr.order,jpeg=TRUE)
	}
}


fam.hm(fam='LTR2B',order.by='CEBPB-ENCFF000YIC')
fam.hm(fam='LTR2C',order.by='TAL1-ENCFF998YDA')
fam.hm(fam='LTR5B',order.by='SPI1-ENCFF000QED')
fam.hm(fam='LTR5_Hs',order.by='SPI1-ENCFF000QED')
fam.hm(fam='LTR12C',order.by='PKNOX1-ENCFF921VOK')
fam.hm(fam='LTR13A',order.by='PKNOX1-ENCFF921VOK')


