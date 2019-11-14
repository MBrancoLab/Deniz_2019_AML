setwd('~/Desktop/AML/AML_scripts/K562_ENCODE')


##ENCODE metadata for all K562 ChIP-seq experiments

meta = read.delim('K562_metadata.tsv',as.is=T)


##get peak bed files

bed = meta[grepl('^bed',meta$File.format) & meta$Output.type=='peaks',]

hg19 = bed[bed$Assembly=='hg19',]
hg38 = bed[bed$Assembly=='GRCh38',]

hg19.targets = unique(hg19$Experiment.target)
hg38.targets = unique(hg38$Experiment.target)

sum(!(hg38.targets %in% hg19.targets))  #2 targets are not in hg19
sum(!(hg19.targets %in% hg38.targets))  #42 targets are not in hg38
#I will get all hg19 and overlap with hg19 annotation of LTRs


##write

df = data.frame(target=gsub('-human','',hg19$Experiment.target),url=hg19$File.download.URL)
write.table(df,'K562_hg19_peaks.txt',sep='\t',quote=F,row.names=F,col.names=F)
