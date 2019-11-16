##Plots average gene expression values for each condition
##Highlights genes that are dCas9 on- or off-targets


setwd('~/Deniz_2019_AML/CRISPRi')


##rna-seq data

rna = read.delim('LTR2B_RNA_vsd.txt',as.is=T)
K562.ctrl = rowMeans(rna[,2:3])
K562.krab = rowMeans(rna[,4:6])
OCI.ctrl = rowMeans(rna[,7:8])
OCI.krab = rowMeans(rna[,9:11])

expr = data.frame(Gene=rna$X,K562.ctrl,K562.krab,OCI.ctrl,OCI.krab)


##add nearest Cas9 peak annotation

anno = read.delim('genes_closest_dCas9.txt',header=F,as.is=T)

match.id = match(expr$Gene,anno$V4)
expr$peak = anno$V10[match.id]
expr$peak.dist = anno$V18[match.id]
expr$LTR = anno$V17[match.id]


##genes within 50kb of a peak

near.ltr = expr$peak.dist<50000 & grepl('LTR2',expr$LTR)
near.off = expr$peak.dist<50000 & !grepl('LTR2',expr$LTR)


##scatter plots

par(mfrow=c(1,2))

plot(expr$K562.ctrl,expr$K562.krab,pch=19,cex=0.3,col='grey',
 xlab='no gRNAs',ylab='LTR2B gRNAs',las=1,main='K562')
points(expr$K562.ctrl[near.off],expr$K562.krab[near.off],pch=19,cex=0.3,col='black')
points(expr$K562.ctrl[near.ltr],expr$K562.krab[near.ltr],pch=19,cex=0.3,col='orange')
abline(a=0,b=1,lty=2)

plot(expr$OCI.ctrl,expr$OCI.krab,pch=19,cex=0.3,col='grey',
 xlab='no gRNAs',ylab='LTR2B gRNAs',las=1,main='OCI-AML3')
points(expr$OCI.ctrl[near.off],expr$OCI.krab[near.off],pch=19,cex=0.3,col='black')
points(expr$OCI.ctrl[near.ltr],expr$OCI.krab[near.ltr],pch=19,cex=0.3,col='orange')
abline(a=0,b=1,lty=2)


##cross with differential expression info

de.k562 = read.table('LTR2B_de_k562.txt',row.names=1)
de.oci = read.table('LTR2B_de_oci.txt',row.names=1)
de = list(de.k562,de.oci)
names(de) = c('K562','OCI3')

de.ltr = lapply(de,function(x) x[rownames(x) %in% expr$Gene[near.ltr],])
de.off = lapply(de,function(x) x[rownames(x) %in% expr$Gene[near.off],])



