<<<<<<< HEAD

#### from Wei Zhou
#https://github.com/globalbiobankmeta/Loci_Definition/blob/main/annotate_novel_or_known.sh

#extend 500kb around the lead variants in the GWAS summary statistics
#merge overlapped regions
#extend 500kb around the previously reported top variants
#annotate the GWAS loci as novel or known based on the overlaps

summary_stats_file=$1
outprefix=$2
##This script will generate regions ${outprefix}.regions.txt which contains merged regions with "knownRegionFlank" bp up- and down-stream top hits

Rscript /net/snowwhite/home/bwolford/2021_analysis/vte/gbmi_vte/ManhattanPlot.r \
	--input=${summary_stats_file}	\
	--PVAL="inv_var_meta_p"	\
	--knownRegionFlank=500000	\
	--prefix=${outprefix}	\
	--ismanhattanplot=TRUE	\
	--isannovar=FALSE	\
	--isqqplot=FALSE	\
	--CHR="#CHR"	\
	--POS=POS	\
	--ALLELE1=REF	\
	--ALLELE2=ALT

############

#region_file=${outprefix}.regions.txt  #region.txt output by Manhattan plot script
#reported_bed_file=$2 #bed file with three columns that contain the region around the previously reported top hits (500kb up- and down-stream)

#tail -n +2 $region_file awk '{if ($3 >= 0) print "chr"$2"\t"$3"\t"$4; else print "chr"$2"\t0\t"$4;}' | sort -k1,1 -k2,2n > ${region_file}.reformat.bed

#bedtools merge -i ${region_file}.reformat.bed > ${region_file}.reformat.mergeOverlap.bed 

#bedtools intersect -a ${region_file}.reformat.mergeOverlap.bed -b ${reported_bed_file} -wao | awk '{ if ($4 == ".") print $1"\t"$2"\t"$3"\tPotentiallyNovel"; else print $1"\t"$2"\t"$3"\tPreviouslyReported";}' | sed 's/chr//g' | sort -k1,1 -k2,2 -g | uniq| awk '{print $0"\t"NR}' > ${region_file}.known_or_novel.bed
=======
region_file=$1  #region.txt output by define_sig_region.sh
reported_bed_file=$2 #bed file with three columns that contain the region around the previously reported top hits (500kb up- and down-stream)

tail -n +2 $region_file awk '{if ($3 >= 0) print "chr"$2"\t"$3"\t"$4; else print "chr"$2"\t0\t"$4;}' | sort -k1,1 -k2,2n > ${region_file}.reformat.bed

bedtools merge -i ${region_file}.reformat.bed > ${region_file}.reformat.mergeOverlap.bed 

bedtools intersect -a ${region_file}.reformat.mergeOverlap.bed -b ${reported_bed_file} -wao | awk '{ if ($4 == ".") print $1"\t"$2"\t"$3"\tPotentiallyNovel"; else print $1"\t"$2"\t"$3"\tPreviouslyReported";}' | sed 's/chr//g' | sort -k1,1 -k2,2 -g | uniq| awk '{print $0"\t"NR}' > ${region_file}.known_or_novel.bed
>>>>>>> 5a56f9ee71bddd4f8693a4c2ab63b52915d6806b
