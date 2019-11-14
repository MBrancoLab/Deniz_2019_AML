##Look for enriched TF motifs in DHS+ vs. DHS- ERVs
##Uses AME tool from MEME and the HOCOMOCO v11 database
##Input are fasta files containing sequences of DHS+ and DHS- ERVs, separated by family


for te in LTR2B LTR2C LTR5B LTR5_Hs LTR12C LTR13A; do
	ame --control ${te}_nodhs.fa --o ${te}_hocomoco ${te}_dhs.fa ~/Documents/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
done
