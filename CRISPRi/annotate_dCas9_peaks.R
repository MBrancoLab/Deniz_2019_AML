##Annotates dCas9 ChIP-seq peaks with LTRs and genic features
##Last line uses bedtools


setwd('~/Deniz_2019_AML/CRISPRi')


##peak data

cas9.file = 'dCas9_peaks.narrowPeak'


##annotations

anno.files = list(ltr2b='../Annotations/LTR2B_hg38_id.bed',
 ltr2='../Annotations/LTR2_hg38_id.bed',
 prom='../Annotations/promoters_hg38.bed',
 exon='../Annotations/exons_hg38.bed',
 intron='../Annotations/introns_hg38.bed')


##intersect peaks with annotations

path_to_intersectbed = '~/Documents/bedtools-2.26.0/bin/intersectBed'

anno = lapply(anno.files, function(x) {
	out=tempfile()
	
	command=paste(path_to_intersectbed,"-c -a", cas9.file,"-b",x,">",out,sep=" ")
	cat(command,"\n")
	try(system(command))
	res=read.table(out,header=F,as.is=T)

	unlink(out)
	return(res)
})


##number of peaks at each feature (some peaks in more than one feature)

unlist(lapply(anno,function(x) sum(x$V11>0)))


##annotate Cas9 peaks
##uses annotation hierarchy: LTR2B/LTR2>promoter>exon>intron>intragenic

cas9 = read.delim(cas9.file,header=F,as.is=T)
cas9$anno = character(nrow(cas9))

cas9$anno[anno$ltr2b$V11>0] = 'LTR2B'
cas9$anno[anno$ltr2$V11>0] = 'LTR2'  #there are no peaks that overlap both an LTR2B and an LTR2 element
cas9$anno[anno$prom$V11>0 & cas9$anno==''] = 'Promoter'
cas9$anno[anno$exon$V11>0 & cas9$anno==''] = 'Exon'
cas9$anno[anno$intron$V11>0 & cas9$anno==''] = 'Intron'
cas9$anno[cas9$anno==''] = 'Intragenic'

write.table(cas9,'dCas9_peaks_anno.txt',sep='\t',quote=F,col.names=F,row.names=F)


##number of peaks at each feature (each peak only counted once)

n.feat = tapply(cas9$anno,factor(cas9$anno),length)


##plot

barplot(rbind(n.feat[5],n.feat[4],n.feat[6],n.feat[1],n.feat[3],n.feat[2]),
 col=c('lightblue','navy','orange','red','green3',grey(0.8)),
 ylim=c(0,400),names.arg='',ylab='n dCas9 peaks',xlim=c(0,2))



#######################################

##annotate genes with nearest Cas9 peak

path_to_closestbed='~/Documents/bedtools-2.26.0/bin/closestBed'
system(paste(path_to_closestbed,'-d -a ../Annotations/GRCh38p12_genes.bed -b dCas9_peaks_anno.txt > genes_closest_dCas9.txt'))


