setwd('~/Desktop/AML/AML_scripts/ChIP-seq')


##file metadata

meta = read.delim('~/Desktop/AML/AML_scripts/blueprint_files.tsv',as.is=T)


##intersect function

intersectBED <- function(a.file,b.file,opt.string="-c") {
	out=tempfile()
	
	command=paste("~/Documents/bedtools-2.26.0/bin/intersectBed -a", a.file,"-b",b.file,opt.string,">",out,sep=" ")
	cat(command,"\n")
	try(system(command))
	res=read.table(out,header=F,as.is=T)

	unlink(out)
	return(res)
}


##overlap function

get.overlaps = function(peak.files,ltr.file,histone) {
	hist.files = peak.files[grep(histone,peak.files)]
	fnames = unlist(lapply(strsplit(hist.files,split='/'),function(x) x[length(x)]))
	
	overlap = intersectBED(ltr.file,paste(hist.files[1]))
	for (i in 2:length(hist.files)) {
		overlap = cbind(overlap,intersectBED(ltr.file,hist.files[i])$V7)
	}
	
	cell = meta$Sub.group[match(fnames,unlist(lapply(strsplit(meta$URL,'/'),function(x) x[length(x)])))]
	cell[cell=='Acute Myeloid Leukemia'] = 'AML'
	cell[cell=='macrophage'] = 'Macro'
	cell[cell=='CD14-positive, CD16-negative classical monocyte'] = 'Mono'	
	snames = paste(gsub('ERX[[:print:]]+','',fnames),cell,sep='')
	
	colnames(overlap) = c('chr','start','end','name','null','strand',snames)
	
	write.table(overlap,paste('./LTR_histone_overlaps/allLTR_',histone,'_overlaps.txt',sep=''),
	 sep='\t',quote=FALSE,row.names=FALSE)
}


##get overlaps for each histone mark

peak.files = list.files('./peak_files',full.names=T)
ltr.file = '../Annotations/all_LTRs.bed'

get.overlaps(peak.files,ltr.file,'H3K27ac')
get.overlaps(peak.files,ltr.file,'H3K4me1')
get.overlaps(peak.files,ltr.file,'H3K4me3')
get.overlaps(peak.files,ltr.file,'H3K9me3')

