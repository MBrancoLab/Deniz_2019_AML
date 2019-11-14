setwd('~/Desktop/AML/AML_scripts/TF_motifs')


##read in ame files

shuf.files = list.files('./ame_files',pattern='shuf_ame.txt')
dhs.files = list.files('./ame_files',pattern='dhs_ame.txt')

shuf = list()
dhs = list()

for (i in 1:length(shuf.files)) {
	shuf[[i]] = scan(paste('./ame_files',shuf.files[i],sep='/'),character(),sep='\n')
	dhs[[i]] = scan(paste('./ame_files',dhs.files[i],sep='/'),character(),sep='\n')
}
names(shuf) = gsub('_shuf_ame.txt','',shuf.files)
names(dhs) = gsub('_dhs_ame.txt','',dhs.files)


##get motifs

get.motifs = function(x) {
	if (length(x)<=8) {
		return()
	} else {
		sp = strsplit(x[9:length(x)],split=' ')
		motifs = unlist(lapply(sp,function(y) y[8]))
		return(motifs)
	}
}

shuf.motifs = lapply(shuf,get.motifs)
dhs.motifs = lapply(dhs,get.motifs)

motifs = sort(unique(unlist(c(shuf.motifs,dhs.motifs)))) ##370 motifs


##make summary table

is.sig = matrix(nrow=length(motifs),ncol=length(shuf.motifs)*2)

for (i in 1:length(shuf.motifs)) {
	is.sig[,i] = motifs %in% shuf.motifs[[i]]
	is.sig[,i+length(shuf.motifs)] = motifs %in% dhs.motifs[[i]]
}
rownames(is.sig) = motifs
colnames(is.sig) = c(paste(names(shuf.motifs),'shuf',sep='_'),paste(names(dhs.motifs),'dhs',sep='_'))

write.table(is.sig,'enriched_motifs.txt',sep='\t',quote=F,col.names=NA)