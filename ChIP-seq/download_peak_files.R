##Downloads relevant histone ChIP-seq peak data from Blueprint

setwd('~/Deniz_2019_AML/ChIP-seq')

dlist = read.delim('../blueprint_files.tsv',as.is=T)


##select cell types of interest

cell = dlist[dlist$Sub.group=='Acute Myeloid Leukemia'|
 dlist$Sub.group=='CD14-positive, CD16-negative classical monocyte' |
 dlist$Sub.group=='macrophage',]


##select histone marks of interest

k27ac = cell[cell$Format=='BED' & cell$Experiment=='H3K27ac',]
k4me1 = cell[cell$Format=='BED' & cell$Experiment=='H3K4me1',]
k4me3 = cell[cell$Format=='BED' & cell$Experiment=='H3K4me3',]
k9me3 = cell[cell$Format=='BED' & cell$Experiment=='H3K9me3',]


##only keep samples profiled for all marks

k27ac = k27ac[k27ac$Donor %in% k4me1$Donor &
 k27ac$Donor %in% k4me3$Donor &
 k27ac$Donor %in% k9me3$Donor,]

k4me1 = k4me1[k4me1$Donor %in% k27ac$Donor &
 k4me1$Donor %in% k4me3$Donor &
 k4me1$Donor %in% k9me3$Donor,]

k4me3 = k4me3[k4me3$Donor %in% k27ac$Donor &
 k4me3$Donor %in% k4me1$Donor &
 k4me3$Donor %in% k9me3$Donor,]

k9me3 = k9me3[k9me3$Donor %in% k27ac$Donor &
 k9me3$Donor %in% k4me1$Donor &
 k9me3$Donor %in% k4me3$Donor,]

hist = list(k27ac,k4me1,k4me3,k9me3)


##remove duplicate samples (keep first one)

dedup = lapply(hist, function(x) {
	donor = factor(x$Donor)
	return(tapply(x$URL,donor,function(y) y[1]))
})


##Get data

for (file in unlist(dedup)) {
	try(system(paste('wget',file)))
}


