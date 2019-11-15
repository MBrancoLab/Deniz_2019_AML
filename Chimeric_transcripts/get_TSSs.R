##Extracts coordinates of TSSs from GTF files

##list GTF files

gtf.files = list.files(pattern='_multi.gtf') ##from remove_single_exon.R


##get TSSs

tss = list()
for (i in 1:length(gtf.files)) {
	gtf = read.delim(gtf.files[i],header=F,as.is=T,quote="")
	gtf = rbind(gtf,data.frame(V1='chrA',V2='MadeUp',V3='exon',V4=1,V5=100,V6=1000,V7='+',V8='.',V9='exon_number "1"',stringsAsFactors=F))
	#added extra bogus line to make it easier to deal with last gene

	exon1 = grep('exon_number "1"',gtf$V9) #5'-most in DNA, not RNA (i.e., mind minus strand genes)
	exon1.rev = which(gtf$V7[exon1]=='-')	
	exon1[exon1.rev] = exon1[exon1.rev+1]-1 #exon1 of rev genes becomes the one before exon1 of the next gene
	
	pos = gtf$V4[exon1]
	pos[exon1.rev] = gtf$V5[exon1[exon1.rev]]
	tss[[i]] = data.frame(chr=gtf$V1[exon1],pos)
}
names(tss) = gsub('_multi.gtf','',gtf.files)

save(tss,file='Blueprint_TSSs.Rdata')
