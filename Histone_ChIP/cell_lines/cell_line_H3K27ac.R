##Finds overlaps between LTRs of interest and H3K27ac peaks from cell lines
##Also processes overlaps with K562 ChromHMM ENCODE data
##Checks the above overlaps against H3K27ac-marked LTRs in Blueprint AML samples
##Uses bedtools

setwd('~/Deniz_2019_AML/ChIP-seq/cell_lines')
path_to_intersectbed = '~/Documents/bedtools-2.26.0/bin/intersectBed'


##### Cell line Data ######


##intersect function

intersectBED <- function(a.file,b.file,opt.string="-u") {
	out=tempfile()
	
	command=paste(path_to_intersectbed,"-a", a.file,"-b",b.file,opt.string,">",out,sep=" ")
	cat(command,"\n")
	try(system(command))
	res=read.table(out,header=F,as.is=T)

	unlink(out)
	return(res)
}


##get overlap between H3K27ac data and LTRs

peak.files = list.files(pattern='broadPeak')
ltr.file = '../../Annotations/all_LTRs.bed'

k27 = list()
for (i in 1:length(peak.files)) {
	k27[[i]] = intersectBED(ltr.file,peak.files[i])
}
names(k27) = gsub('\\.broadPeak','',peak.files)



##### Blueprint AML Data ######


bp = read.delim('../LTR_histone_overlaps/allLTR_H3K27ac_overlaps.txt',as.is=T)

aml.sum = rowSums(bp[,grep('AML',colnames(bp))])

in.bp = lapply(k27,function(x) {
	aml0 = sum(x$V4 %in% bp$name[aml.sum==0])/sum(aml.sum==0)
	aml5 = sum(x$V4 %in% bp$name[aml.sum>=5 & aml.sum<10])/sum(aml.sum>=5 & aml.sum<10)
	aml10 = sum(x$V4 %in% bp$name[aml.sum>10])/sum(aml.sum>10)
	return(c(aml0,aml5,aml10))
})

in.bp = matrix(unlist(in.bp),ncol=length(in.bp))
colnames(in.bp) = names(k27)



##### HMM Data ######
#From ENCODE, hg38 ChromHMM file in Annotations folder

hmm = read.delim('ChromHMM_LTR_overlaps.txt',header=F,as.is=T) #overlaps generated with intersectBed


##calculate % of LTR overlap

p.over = function(x) {
	bp = x$V9-x$V8
	left = x$V8-x$V2
	right = x$V3 - x$V9
	if (left>=0 & right>=0) {
		over = 1
	} else {
		bp.over = bp
		if (left<0) bp.over = bp.over+left
		if (right<0) bp.over = bp.over+right
		over = bp.over/bp
	}
	return(over)
}

ltr.over = numeric(nrow(hmm))
for (i in 1:nrow(hmm)) ltr.over[i] = p.over(hmm[i,])


##remove duplicate LTR entries based on % overlap

hmm2 = hmm[order(hmm$V10,ltr.over),]

to.remove = logical(nrow(hmm2))
for (i in 2:nrow(hmm)) {
	if (hmm2$V10[i]==hmm2$V10[i-1]) {
		to.remove[i-1] = TRUE
	}
}

hmm.dedup = hmm2[!to.remove,]


##count enhancers/heterochromatin in AML LTRs

is.enh = grepl('Enhancer',hmm.dedup$V4)
is.pro = grepl('Promoter',hmm.dedup$V4)

enh = c(sum(is.enh & hmm.dedup$V10 %in% bp$name[aml.sum==0])/sum(aml.sum==0),
 sum(is.enh & hmm.dedup$V10 %in% bp$name[aml.sum>=5 & aml.sum<10])/sum(aml.sum>=5 & aml.sum<10),
 sum(is.enh & hmm.dedup$V10 %in% bp$name[aml.sum>10])/sum(aml.sum>10))

pro = c(sum(is.pro & hmm.dedup$V10 %in% bp$name[aml.sum==0])/sum(aml.sum==0),
 sum(is.pro & hmm.dedup$V10 %in% bp$name[aml.sum>=5 & aml.sum<10])/sum(aml.sum>=5 & aml.sum<10),
 sum(is.pro & hmm.dedup$V10 %in% bp$name[aml.sum>10])/sum(aml.sum>10))



##### Plot ######


barplot(cbind(in.bp,enh,pro)*100,beside=T,las=2,cex.names=0.5,ylim=c(0,80),
 col=c('azure','lightblue','navy'),ylab='% of elements')


