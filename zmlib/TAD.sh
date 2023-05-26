#----------------------------------------------
#ZO1
SAM="WT"
mkdir -p ${path}/TAD/Dense
mkdir -p ${path}/TAD/DI_mat
mkdir -p ${path}/TAD/DI
python ${Script}/sparseToDense.py -b ${Script}/abs/hg19_20000_abs.bed /media/wulab/C106/2018-7-10_HiC/HiC-Pro/hic_results/matrix/${SAM}/iced/20000/${SAM}_20000_iced.matrix --perchr -o ${path}/TAD/Dense/sample

chroms=`cut -f 1 ${Script}/abs/hg19_10000_abs.bed|uniq |grep -v "chrM"`
for chr in $chroms;do
	echo $chr
	python ${Script}/DI/hic_di.py -a ${Script}/abs/hg19_20000_abs.bed -m ${path}/TAD/Dense/sample_${chr} -c ${chr} -o ${path}/TAD/DI_mat/sample_${chr}_di.matrix
	perl ${Script}/DI/DI_from_matrix.pl ${path}/TAD/DI_mat/sample_${chr}_di.matrix 20000 2000000 ${Script}/DI/hg19.chrom.sizes > ${path}/TAD/DI/sample_DI_${chr}.bed
done
cat ${path}/TAD/DI/sample_DI_*.bed > ${path}/TAD/${SAM}_DI_20Kb.bed
cat ${path}/TAD/${SAM}_DI_20Kb.bed|awk -v OFS="\t" '{if($1 == "X") print "23",$2,$3,$4;else if($1 == "Y") print "24",$2,$3,$4; else print $1,$2,$3,$4}' > ${path}/TAD/${SAM}_DI_hmm_20Kb.bed
rm -r ${path}/TAD/Dense && rm -r ${path}/TAD/DI_mat && rm -r ${path}/TAD/DI

nice matlab HMM_calls.m # ${path}/TAD/${SAM}_DI_hmm.bed, ${path}/TAD/${SAM}_DI_hmm_output.bed

perl ${Script}/DI/file_ends_cleaner.pl ${path}/TAD/${SAM}_DI_hmm_20Kb_output.bed ${path}/TAD/${SAM}_DI_hmm_20Kb.bed | perl ${Script}/DI/converter_7col.pl > ${path}/TAD/${SAM}_hmm_7colfile
chroms=`cut -f 1 ${path}/TAD/${SAM}_hmm_7colfile|uniq`
for chr in $chroms;do
	echo $chr
	grep $chr ${path}/TAD/${SAM}_hmm_7colfile > ${path}/TAD/colfile_${chr}.bed
	perl ${Script}/DI/hmm_probablity_correcter.pl ${path}/TAD/colfile_${chr}.bed 2 0.99 20000 | perl ${Script}/DI/hmm-state_caller.pl faifile ${chr} | perl ${Script}/DI/hmm-state_domains.pl > ${path}/TAD/sample_${chr}.dm
	rm ${path}/TAD/colfile_${chr}.bed
done
cat ${path}/TAD/sample_chr*.dm | awk -v OFS="\t" '{if($1 == "chr23") print "chrX",$2,$3;else if($1 == "chr24") print "chrY",$2,$3; else print $1,$2,$3}'|sort -k1,1 -k2,2n > ${path}/TAD/${SAM}_Domain.bed
rm ${path}/TAD/sample_chr*.dm && rm ${path}/TAD/${SAM}_hmm_7colfile && rm ${path}/TAD/${SAM}_DI_hmm_20Kb_output.bed && rm ${path}/TAD/${SAM}_DI_hmm_20Kb.bed


