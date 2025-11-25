#with SRAtoolkit installed

prefetch SRR8147026 --max-size 420000000000
fasterq-dump SRR8147026 --split-files --include-technical

prefetch SRR8147027 --max-size 420000000000
fasterq-dump SRR8147027 --split-files --include-technical

prefetch SRR8147028 --max-size 420000000000
fasterq-dump SRR8147028 --split-files --include-technical

prefetch SRR8147029 --max-size 420000000000
fasterq-dump SRR8147029 --split-files --include-technical


#concatenate reads
cat SRR8147026_1.fastq SRR8147027_1.fastq SRR8147028_1.fastq SRR8147029_1.fastq > Undetermined_S0_L001_R1_001.fastq

cat SRR8147026_2.fastq SRR8147027_2.fastq SRR8147028_2.fastq SRR8147029_2.fastq > Undetermined_S0_L001_R2_001.fastq

cat SRR8147026_3.fastq SRR8147027_3.fastq SRR8147028_3.fastq SRR8147029_3.fastq > Undetermined_S0_L001_R3_001.fastq

cat SRR8147026_4.fastq SRR8147027_4.fastq SRR8147028_4.fastq SRR8147029_4.fastq > Undetermined_S0_L001_R4_001.fastq

pigz -p 4 Undetermined_S0_L001_R1_001.fastq
pigz -p 4 Undetermined_S0_L001_R2_001.fastq
pigz -p 4 Undetermined_S0_L001_R3_001.fastq
pigz -p 4 Undetermined_S0_L001_R4_001.fastq


