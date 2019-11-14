##Overlaps DHS data with Repeatmasker annotation
##Generates 1000 shuffled versions of the DHS data within mappable regions of the genome
##Compares the real overlap with the random overlaps to identify significantly enrinched/depleted repeat families

##Note: DHS file name is passed on as an argument from command line

##required files

source('rmsk_intersect_functions.R')  ##make.random, intersectBED and test.overlap functions
rmasker = 'Repeatmasker_hg38.txt' ##Repeatmasker annotation file in bed format, uncompressed (not in this repository)
rmasker.counts = 'hg38_rm_counts.txt' ##Total counts of each repeatmasker family (repName)
unmappable = 'unmappable_hg38_MB.bed'  ##Unmappable regions of the genome in bed format, uncompressed


##read name of bed file containing DHS data

args = commandArgs(trailingOnly=TRUE)
npf.file = args[1]


##intersect with repeatmasker

npf = read.delim(npf.file,header=F,as.is=T)
rm.count = read.table(rmasker.counts)

bed = npf[npf$V3-npf$V2<5000,1:3]  ##exclude enormous peaks
bed.rm = intersectBED(bed,rm.file=rmasker)
bed.rm = unique(bed.rm)
real.overlap = tapply(bed.rm$V5,factor(bed.rm$V5,levels=levels(rm.count$V1)),length)
real.overlap[is.na(real.overlap)] = 0

n=1000
rand.overlaps = matrix(nrow=length(real.overlap),ncol=n)
for (i in 1:n) {
	bed.random = make.random(bed,unmap.file=unmappable,genome='hg38')
	bed.rand.rm = intersectBED(bed.random,rm.file=rmasker)
	bed.rand.rm = unique(bed.rand.rm)
	rand.overlaps[,i] = tapply(bed.rand.rm$V5,factor(bed.rand.rm$V5,levels=levels(rm.count$V1)),length)
	rand.overlaps[is.na(rand.overlaps[,i]),i] = 0
}


##calculate p value

left.tail = rowSums(rand.overlaps<=as.vector(real.overlap))/n
right.tail = rowSums(rand.overlaps>=as.vector(real.overlap))/n

tails = cbind(left.tail,right.tail)
pval = apply(tails,1,function(x) min(x)*2)
pval[pval>1] = 1


##get enrichment

random.mean = rowMeans(rand.overlaps)
log.enrich = log2(real.overlap/random.mean)
total = rm.count$V2[match(names(real.overlap),rm.count$V1)]

out = data.frame(log.enrich,pval,total,real.overlap,random.mean)
write.table(out,paste(npf.file,'overlap.txt',sep='_'),sep='\t',quote=F,col.names=NA)
