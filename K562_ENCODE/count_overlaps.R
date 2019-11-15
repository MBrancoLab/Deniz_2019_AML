##Counts total number of LTRs from each family that overlaps TF ChIP-seq peak data from K562
##Generates a count matrix
##Uses bedtools


setwd('~/Deniz_2019_AML/K562_ENCODE')
path_to_intersectbed = '~/Documents/bedtools-2.26.0/bin/intersectBed'


##intersectBED function

intersectBED <- function(a.file,b.file,opt.string="-u") {
	out=tempfile()
	
	command=paste(path_to_intersectBed,"-a", a.file,"-b",b.file,opt.string,">",out,sep=" ")
	cat(command,"\n")
	try(system(command))
	res=read.table(out,header=F,as.is=T,col.names=paste('V',1:6,sep=''))

	unlink(out)
	return(res)
}


##TF and LTR files

bed.files = list.files('./peak_files',full.names=TRUE) #not included in the repository
ltr.files = list.files('./LTRs_hg19',full.names=TRUE) #change to run through shuffled LTR files


##get overlaps

count = matrix(nrow=length(bed.files),ncol=length(ltr.files))
for (i in 1:length(bed.files)) {
	for (j in 1:length(ltr.files)) {
		over = intersectBED(ltr.files[j],bed.files[i])
		count[i,j] = nrow(over)
	}
}
rownames(count) = bed.files
colnames(count) = gsub('_hg19.bed','',ltr.files)

write.table(count,'LTR_hg19_overlaps.txt',sep='\t',quote=F,col.names=NA)