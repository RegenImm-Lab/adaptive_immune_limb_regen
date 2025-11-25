#with SRAtoolkit installed

prefetch SRR8147030 --max-size 420000000000
fasterq-dump SRR8147030 --split-files --include-technical

prefetch SRR8147031 --max-size 420000000000
fasterq-dump SRR8147031 --split-files --include-technical

prefetch SRR8147032 --max-size 420000000000
fasterq-dump SRR8147032 --split-files --include-technical

prefetch SRR8147033 --max-size 420000000000
fasterq-dump SRR8147033 --split-files --include-technical

#concatenate reads
cat SRR8147030_1.fastq SRR8147031_1.fastq SRR8147032_1.fastq SRR8147033_1.fastq > Undetermined_S0_L001_R1_001.fastq

cat SRR8147030_2.fastq SRR8147031_2.fastq SRR8147032_2.fastq SRR8147033_2.fastq > Undetermined_S0_L001_R2_001.fastq

cat SRR8147030_3.fastq SRR8147031_3.fastq SRR8147032_3.fastq SRR8147033_3.fastq > Undetermined_S0_L001_R3_001.fastq

cat SRR8147030_4.fastq SRR8147031_4.fastq SRR8147032_4.fastq SRR8147033_4.fastq > Undetermined_S0_L001_R4_001.fastq

pigz -p 4 Undetermined_S0_L001_R1_001.fastq
pigz -p 4 Undetermined_S0_L001_R2_001.fastq
pigz -p 4 Undetermined_S0_L001_R3_001.fastq
pigz -p 4 Undetermined_S0_L001_R4_001.fastq


