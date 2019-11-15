##Extracts a list of bigwig files that correspond to the histone ChIP-seq files used in the analysis

setwd('~/Deniz_2019_AML/')


##blueprint files

dlist = read.delim('blueprint_files.tsv',as.is=T)


##peak files

peak.files = list.files(path='./ChIP-seq/peak_files')


##select bigwig files for the same experiments

serx = gsub('\\.bwa.GRCh38[[:print:]]+','',peak.files)

bw = character(length(serx))
for (i in 1:length(serx)) {
	bw[i] = dlist$URL[grepl(serx[i],dlist$URL) & dlist$Format=='bigWig']
}


##write

write(bw,'bw_files.txt')
