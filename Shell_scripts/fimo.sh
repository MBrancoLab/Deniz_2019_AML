##Identifies and annotates TF motifs in individual ERVs
##Uses FIMO tool from MEME and the HOCOMOCO v11 database
##Input are fasta files separated by family

database=~/Documents/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme

for te in LTR2B LTR2C LTR5B LTR5_Hs LTR12C LTR13A; do
	fimo --o ${te}_fimo $database ../TF_motifs/fasta_files/${te}_dhs.fa 
done
