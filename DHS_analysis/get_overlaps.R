##Generates a matrix of DHS-TE overlaps for the selected TE families
##Uses bedtools

setwd('~/Deniz_2019_AML/DHS_analysis')
path_to_intersectbed = "~/Documents/bedtools-2.26.0/bin/intersectBed"


##intersect function

intersectBED <- function(a.file,b.file,opt.string="-c") {
	out=tempfile()
	
	command=paste(path_to_intersectbed,"-a", a.file,"-b",b.file,opt.string,">",out,sep=" ")
	cat(command,"\n")
	try(system(command))
	res=read.table(out,header=F,as.is=T)

	unlink(out)
	return(res)
}


##get peak and ltr data

dhs.files = list.files('./peak_files',full.names=T)
dhs.files = dhs.files[!grepl('Macrophages.C005VG45',dhs.files)] #remove macrophage sample with very few peaks

fnames = unlist(lapply(strsplit(dhs.files,split='/'),function(x) x[length(x)]))

ltr.file = '../Annotations/all_LTRs.bed'  #Repeatmasker annotations for LTR2B, LTR2C, LTR5B, LTR5_Hs, LTR12C, LTR13A


##overlap

overlap = intersectBED(ltr.file,paste(dhs.files[1]))
for (i in 2:length(dhs.files)) {
	overlap = cbind(overlap,intersectBED(ltr.file,dhs.files[i])$V7)
}

snames = gsub('-hg38-unique.npf','',fnames)
snames = gsub('_DHS-unique.npf','',snames)
snames = gsub('_filt.npf','',snames)
snames = gsub('\\.ERX[[:print:]]+','',snames)

colnames(overlap) = c('chr','start','end','name','null','strand',snames)

write.table(overlap,'allLTR_DHS_overlaps.txt',sep='\t',quote=FALSE,row.names=FALSE)
