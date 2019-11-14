setwd('~/Desktop/AML/AML_scripts/TF_motifs/fasta_files')


##DHS data

dhs = read.delim('~/Desktop/AML/AML_scripts/DHS_analysis/allLTR_DHS_overlaps.txt',as.is=T)


##select elements with DHS in at least 5 samples

dhs5 = dhs$name[rowSums(dhs[,7:ncol(dhs)])>=5]


##function to split fasta file

split.fa = function(fam,dhs.ltrs) {
	fa = scan(paste(fam,'hg38.fa',sep='_'),character(),sep='\n')
	new.seq = grep('>',fa)
	for (i in 1:length(new.seq)) {
		ltr = gsub('>','',fa[new.seq[i]])
		seq = fa[new.seq[i]:(new.seq[i]+1)]
		if (ltr %in% dhs.ltrs) {
			write(seq,paste(fam,'dhs.fa',sep='_'),append=TRUE)
		} else {
			write(seq,paste(fam,'nodhs.fa',sep='_'),append=TRUE)
		}
	}
}


##split files

split.fa('LTR2B',dhs5)
split.fa('LTR2C',dhs5)
split.fa('LTR5B',dhs5)
split.fa('LTR5_Hs',dhs5)
split.fa('LTR12C',dhs5)
split.fa('LTR13A',dhs5)