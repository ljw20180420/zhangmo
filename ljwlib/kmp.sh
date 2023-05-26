#!/bin/bash
chr=$1
start=$2
end=$3
excu='LPM/build/LPM'
files='K2_rep1_R1.fq.gz K2_rep1_R2.fq.gz K2_rep2_R1.fq.gz K2_rep2_R2.fq.gz K45_rep1_R1.fq.gz K45_rep1_R2.fq.gz K45_rep2_R1.fq.gz K45_rep2_R2.fq.gz AID_ZO_rep1_R1.fq.gz AID_ZO_rep1_R2.fq.gz AID_ZO_rep2_R1.fq.gz AID_ZO_rep2_R2.fq.gz'
threshold=20
# for file in $files; do
# 	(zcat hic_rawdata/$file | awk -v OFS="\t" '{if(NR%4==2) print 0, $0;}' | $excu znf143.fa forward $threshold > check_genome_run/${file}.LPM.for1; awk -v OFS="\t" '{print $1, $NF;}' check_genome_run/${file}.LPM.for1 | $excu znf143.fa forward 0 > check_genome_run/${file}.LPM.for2) &
# done
# wait

for file in $files; do
	awk -v OFS="\t" -v chr=$chr -v start=$start -v for2=check_genome_run/${file}.LPM.for2 -v threshold=$threshold '
	BEGIN{
		PAIRNUM=0;
	}
	{
		getline line < for2;
		n1 = split($0, fields1, "\t");
		n2 = split(line, fields2, "\t");
		if (n2==2 || fields2[1]<threshold)
			next;
		score = 1 / ((n1-2)/2 * (n2-2)/2);
		for (i=2; i<n1; i+=2)
			for (j=2; j<n2; j+=2)
			{
				start1 = start+fields1[i]
				end1 = start+fields1[i]+fields1[1]
				start2 = start+fields2[j]
				end2 = start+fields2[j]+fields2[1]
				if(fields1[i+1]=="-")
				{
					tmp = start1
					start1 = end1
					end1 = tmp
				}
				if(fields2[j+1]=="-")
				{
					tmp = start2
					start2 = end2
					end2 = tmp
				}
				print chr, start1, end1, chr, start2, end2, "pair" ++PAIRNUM, score, fields1[i+1], fields2[j+1];
			}
	}
	' check_genome_run/${file}.LPM.for1 > check_genome_run/${file}.LPM.for &
done
wait