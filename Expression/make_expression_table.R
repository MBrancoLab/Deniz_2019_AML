setwd('~/Desktop/AML/AML_scripts/Expression')


##Blueprint files

dlist = read.delim('../blueprint_files.tsv')


##select RNA-seq data of interest

sel = dlist[dlist$Sub.group=='Acute Myeloid Leukemia'|
 dlist$Sub.group=='CD14-positive, CD16-negative classical monocyte' |
 dlist$Sub.group=='macrophage' |
 dlist$Sub.group=='common myeloid progenitor' |
 dlist$Sub.group=='hematopoietic stem cell',]
rna = sel[sel$File.type=='Transcription quantification (Genes)',]


##download data

for (file in rna$URL) {
	try(system(paste('wget',file)))
}

rna.file = unlist(lapply(strsplit(as.character(rna$URL),'/'),function(x) x[length(x)]))


##merge all expression files

rna.df = read.delim(rna.file[1],as.is=T)[c(1,7)]
for (i in 2:length(rna.file)) {
	rna.df = cbind(rna.df,read.delim(rna.file[i],as.is=T)$FPKM)
}
colnames(rna.df)[-1] = paste(as.character(rna$Sub.group),as.character(rna$Donor),sep='.')


##find nearest LTR

system('~/Documents/bedtools-2.26.0/bin/sortBed -i ../Annotations/all_LTRs.bed > sorted_LTRs.bed')
system('~/Documents/bedtools-2.26.0/bin/closestBed -d -a ../Annotations/GRCh38p12_genes.bed -b sorted_LTRs.bed > closest_LTR.bed')

close = read.delim('closest_LTR.bed',as.is=T,header=F)
unlink(c('sorted_LTRs.bed','closest_LTR.bed'))


##match gene IDs to add extra info
 
stable.id = gsub('\\.[[:digit:]]+','',rna.df$gene_id)
mid = match(stable.id,close$V5)

rna.df$gene_name = close$V4[mid]
rna.df$LTR_name = close$V10[mid]
rna.df$LTR_distance = close$V13[mid]


##write expression table

write.table(rna.df,'BP_RNA_FPKM.txt',sep='\t',quote=F,row.names=F)


