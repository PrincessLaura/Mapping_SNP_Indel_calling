#!/bin/bash

echo "Enter zipped fastq file 1 name and location, followed by [ENTER]:"
read fastq1
echo "Enter zipped fastq file 2 name and location, followed by [ENTER]:"
read fastq2
echo "Enter name and location of reference genome, followed by [ENTER]:"
read reference_genome_fasta
#echo "Enter name and location of reference genome without .fasta, followed by [ENTER]:"
#read reference_genome_fasta2
echo "Enter analysis name, followed by [ENTER]:"
read NAME

gunzip -c $fastq1 | sed 's/ /_/' > $NAME"_Processed_fastq1"
gunzip -c $fastq2 | sed 's/ /_/' > $NAME"_Processed_fastq2"
#cat $NAME"_Processed_fastq1" $NAME"_Processed_fastq2" >> $NAME"_Processed_fastq"
#cat $fastq1 $fastq2 >> $NAME"_Processed_fastq"

java -jar /pub6/kate/lamotrigine/Trimmed/virusfinder/VirusFinder2.0/bin/CreateSequenceDictionary.jar R= $reference_genome_fasta O= $reference_genome_fasta2".dict"

/pub9/laura/bwa-0.7.10/bwa index -a bwtsw $reference_genome_fasta
#/pub9/laura/bwa-0.7.10/bwa index -a is $reference_genome_fasta

/pub9/laura/bwa-0.7.10/bwa mem -t 20 -R '@RG\tID:Unknwn\tPL:Illumina\tLB:library\tSM:Unknown' -M $reference_genome_fasta $fastq1 $fastq2 | awk '$1 ~ /^@/ || $2 == 65 || $2 == 129 || $2 == 67 || $2 == 131 || $2 == 113 || $2 == 177 || $2 == 81 || $2 == 161 || $2 == 163 || $2 == 83 || $2 == 97 || $2 == 145 || $2 == 99 || $2 == 147 || $2 == 137 || $2 == 73 {print $0}'  > $NAME"_filter.sam"
samtools view -u -q 10 -T $reference_genome_fasta $NAME"_filter.sam" | samtools sort - $NAME"_sort"
samtools index $NAME"_sort.bam"

java -jar -Xmx10g /pub15/xliu/software/picard-tools-1.85/MarkDuplicates.jar I= $NAME"_sort.bam" O= $NAME"_remove_dups.bam" M=duplication.txt REMOVE_DUPLICATES=true AS=true

samtools sort $NAME"_remove_dups.bam" $NAME"_remove_dups_sort"

samtools index $NAME"_remove_dups_sort.bam"

perl /pub9/laura/coverageStatsSplitByChr_v2.pl -i $NAME"_remove_dups_sort.bam" > $NAME"_coverage"

awk '{sum=sum+$4} END {print "Average % coverage of reference contigs=" sum/NR}' $NAME"_coverage" > $NAME"_Coverage_stats"
awk '{sum=sum+$5} END {print "Average depth of coverage of reference contigs=" sum/NR}' $NAME"_coverage" >> $NAME"_Coverage_stats"
echo "Number of mapped contigs:" >> $NAME"_Coverage_stats"
wc -l $NAME"_coverage" >> $NAME"_Coverage_stats"

samtools mpileup -f $reference_genome_fasta $NAME"_remove_dups_sort.bam" > $NAME"_raw.pileup"
awk '$4 != 0 {print $0}' $NAME"_raw.pileup" > $NAME"_final.pileup"
echo "Number of mapped bases:" >> $NAME"_Coverage_stats"
wc -l $NAME"_final.pileup" >> $NAME"_Coverage_stats"
echo "Number of mapped bases at 5X or more:" >> $NAME"_Coverage_stats"
awk '$4 >= 5 {print $0}' $NAME"_final.pileup" | wc -l >> $NAME"_Coverage_stats"
echo "Number of mapped bases at 10X or more:" >> $NAME"_Coverage_stats"
awk '$4 >= 10 {print $0}' $NAME"_final.pileup" | wc -l >> $NAME"_Coverage_stats"
echo "Number of mapped reads 1st pass:" >> $NAME"_Coverage_stats"
grep -v "^@" $NAME"_filter.sam" | wc -l >> $NAME"_Coverage_stats"
echo "Number of uniquely mapped reads:" >> $NAME"_Coverage_stats"
samtools view $NAME"_sort.bam" | grep -v "^@" | wc -l >> $NAME"_Coverage_stats"
echo "Number of uniquely mapped reads after remove duplicates:" >> $NAME"_Coverage_stats"
samtools view $NAME"_remove_dups_sort.bam" | grep -v "^@" | wc -l >> $NAME"_Coverage_stats"

java -jar /pub15/xliu/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R $reference_genome_fasta -I $NAME"_remove_dups_sort.bam" -T RealignerTargetCreator -o $NAME"_target.intervals"

java -jar /pub15/xliu/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R $reference_genome_fasta -I $NAME"_remove_dups_sort.bam" -T IndelRealigner -targetIntervals $NAME"_target.intervals" -o $NAME"_cleaned.bam"

samtools sort $NAME"_cleaned.bam" $NAME"_cleaned_sort"

samtools index $NAME"_cleaned_sort.bam"

java -jar /pub15/xliu/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R $reference_genome_fasta -I $NAME"_cleaned_sort.bam" -T UnifiedGenotyper -stand_call_conf 20.00 -o $NAME".raw.vcf"

java -jar /pub15/xliu/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R $reference_genome_fasta  -T VariantFiltration --variant $NAME".raw.vcf" -o $NAME"_vfallelecalls.vcf" --clusterSize 3 --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 || QUAL < 30.0 || QD < 1.5" --filterName "DodgySNPs" --genotypeFilterExpression "isHet == 1" --genotypeFilterName Heterozygote --genotypeFilterExpression "isHomVar == 1" --genotypeFilterName Homozygote --genotypeFilterExpression "isHomRef == 1" --genotypeFilterName HomozygoteRef

perl /pub9/laura/Correct_for_chr3_split.pl $NAME"_vfallelecalls.vcf" > $NAME"_vfallelecalls_corrected.vcf"

java -jar /pub15/xliu/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R $reference_genome_fasta -I $NAME"_cleaned_sort.bam" -T UnifiedGenotyper -stand_call_conf 20.00 -glm INDEL -minIndelCnt 5 -o $NAME".raw_indel.vcf"

java -jar /pub15/xliu/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R $reference_genome_fasta  -T VariantFiltration --variant $NAME".raw_indel.vcf" -o $NAME"_vfindelcalls.vcf" --clusterSize 3 --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 || QUAL < 30.0 || QD < 1.5" --filterName "Dodgyindels"

perl /pub9/laura/Correct_for_chr3_split.pl $NAME"_vfindelcalls.vcf" > $NAME"_vfindelcalls_corrected.vcf"

echo "Final number of mapped reads:" >> $NAME"_Coverage_stats"
samtools view $NAME"_cleaned_sort.bam" | grep -v "^@" | wc -l >> $NAME"_Coverage_stats"
echo "Number of SNPs:" >> $NAME"_Coverage_stats"
grep "PASS" $NAME"_vfallelecalls_corrected.vcf" | wc -l >> $NAME"_Coverage_stats"
echo "Number of homozygous SNPs:" >> $NAME"_Coverage_stats"
grep "PASS" $NAME"_vfallelecalls_corrected.vcf" | grep "Homozygote:" | wc -l >> $NAME"_Coverage_stats"
echo "Number of Indels:" >> $NAME"_Coverage_stats"
grep "PASS" $NAME"_vfindelcalls_corrected.vcf" | wc -l >> $NAME"_Coverage_stats"
