
##Make random control probes
##unmap.file is a bed file of unmappable regions of the genome
##To shuffle fragments anywhere in the genome, make unmap.file=''

make.random = function(bed,unmap.file='hg38_unmappable.bed',genome='hg38',seed=NA) {
	n = nrow(bed)
	roi.length = bed$V3 - bed$V2
	
	if (grepl('hg',genome)) chr.name = paste('chr',c(as.character(1:22),'X','Y'),sep='')
	if (grepl('mm',genome)) chr.name = paste('chr',c(as.character(1:19),'X','Y'),sep='')
	
	chr.length = list(hg19=c(249250621,243199373,198022430,191154276,180915260,171115067,
	 159138663,146364022,141213431,135534747,135006516,133851895,
	 115169878,107349540,102531392,90354753,81195210,78077248,59128983,
	 63025520,48129895,51304566,155270560,59373566),
	  hg38=c(248956422,242193529,198295559,190214555,181538259,170805979,
	 159345973,145138636,138394717,133797422,135086622,133275309,
	 114364328,107043718,101991189,90338345,83257441,80373285,58617616,
	 64444167,46709983,50818468,156040895,57227415),
	  mm9=c(197195432,181748087,159599783,155630120,152537259,149517037,152524553,
	 131738871,124076172,129993255,121843856,121257530,120284312,125194864,
	 103494974,98319150,95272651,90772031,61342430,166650296,15902555),
	  mm10=c(195471971,182113224,152532351,156508116,151834684,149736546,145441459,
	 129401213,124595110,130694993,122082543,120129022,120421639,124902244,
	 104043685,98207768,94987271,90702639,61431566,171031299,91744698))
	
	gen.id = which(names(chr.length)==genome)
	gen = cbind(chr.name,chr.length[[gen.id]])
	gen.file = tempfile()
	write.table(gen,file=gen.file,quote=F,sep="\t",col.names=F,row.names=F)
	
	i.file = tempfile()
	write.table(bed,file=i.file,quote=F,sep="\t",col.names=F,row.names=F)
	
	out.file = tempfile()
	
	if (is.na(seed)) {
		seed = round(runif(1,0,1000000))
	}
	
	if (unmap.file=='') {
		command=paste("shuffleBed -i", i.file,"-g",gen.file,"-seed",seed,"-noOverlapping >",out.file,sep=" ")
	} else {
		command=paste("shuffleBed -i",i.file,"-g",gen.file,"-seed",seed,"-excl",unmap.file,"-f 0 -noOverlapping >",out.file,sep=" ")
	}	
	cat(command,"\n")
	try(system(command))
	res=read.table(out.file,header=F,as.is=T)
	
	unlink(gen.file)
	unlink(i.file)
	unlink(out.file)
	return(res)
}



##intersectBED with Repeatmasker

intersectBED <- function(bed,opt.string="-wa",rm.file='Repeatmasker_hg38.txt') {
	b.file=tempfile()
	write.table(bed,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
	out=tempfile()
	
	command=paste("intersectBed -a", rm.file,"-b",b.file,opt.string,">",out,sep=" ")
	cat(command,"\n")
	try(system(command))
	res=read.table(out,header=F,as.is=T)

	unlink(b.file)
	unlink(out)
	return(res)
}


##test for enrichment/depletion

test.overlap = function(overlap1,overlap2,counts.file='hg38_rm_counts.txt') {
	rm.count = read.table(counts.file)
	overlap1$V5 = factor(overlap1$V5, levels=levels(rm.count$V1))
	overlap2$V5 = factor(overlap2$V5, levels=levels(rm.count$V1))
	
	overlap1.count = tapply(overlap1$V1,overlap1$V5,length)
	overlap1.count[is.na(overlap1.count)] = 0
	overlap2.count = tapply(overlap2$V1,overlap2$V5,length)
	overlap2.count[is.na(overlap2.count)] = 0
	
	total.count = rm.count$V2
	counts = cbind(total.count,overlap1.count,overlap2.count)
	p = apply(counts, 1, function(x) {
		mat = cbind(x[2:3],x[1]-x[2:3])
		p = fisher.test(mat)$p.value
		return(p)
	})
	
	padj = p.adjust(p, method='BH')
	df = data.frame(levels(rm.count$V1),overlap1.count,overlap2.count,total.count,padj)
	colnames(df) = c('Repeat','Overlap1','Overlap2','Total','Padj')
	return(df)
}
