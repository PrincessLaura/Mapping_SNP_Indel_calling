# Mapping_SNP_Indel_calling
Mapping and SNP/Indel calling pipeline implementing BWA and GATK. Designed for use with the IWGSC pseudo chromosomal molecules reference sequence (chr 3B split at 480000000bp to allow compatibility with SAMtools reference indexing).

Requires: GATK, BWA, SAMtools, picard tools-MarkDuplicates, coverageStatsSplitByChr_v2.pl, awk, Correct_for_chr3_split.pl
