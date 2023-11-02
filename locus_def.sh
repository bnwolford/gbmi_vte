region_file=$1  #region.txt output by define_sig_region.sh
reported_bed_file=$2 #bed file with three columns that contain the region around the previously reported top hits (500kb up- and down-stream)

tail -n +2 $region_file awk '{if ($3 >= 0) print "chr"$2"\t"$3"\t"$4; else print "chr"$2"\t0\t"$4;}' | sort -k1,1 -k2,2n > ${region_file}.reformat.bed

bedtools merge -i ${region_file}.reformat.bed > ${region_file}.reformat.mergeOverlap.bed 

bedtools intersect -a ${region_file}.reformat.mergeOverlap.bed -b ${reported_bed_file} -wao | awk '{ if ($4 == ".") print $1"\t"$2"\t"$3"\tPotentiallyNovel"; else print $1"\t"$2"\t"$3"\tPreviouslyReported";}' | sed 's/chr//g' | sort -k1,1 -k2,2 -g | uniq| awk '{print $0"\t"NR}' > ${region_file}.known_or_novel.bed
