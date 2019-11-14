##Identifies and annotates TF motifs in individual ERVs
##Uses FIMO tool from MEME and the HOCOMOCO v11 database
##Input are fasta files separated by family


for te in LTR2B LTR2C LTR5B LTR5_Hs LTR12C LTR13A; do
	fimo --o ${te}_fimo ~/Documents/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme ${te}_dhs.fa 
done
