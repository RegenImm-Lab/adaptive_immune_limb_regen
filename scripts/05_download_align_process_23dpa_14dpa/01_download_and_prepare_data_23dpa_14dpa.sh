#with SRAtoolkit installed
prefetch SRR8147022 --max-size 420000000000
fasterq-dump SRR8147022 --split-files --include-technical

prefetch SRR8147023 --max-size 420000000000
fasterq-dump SRR8147023 --split-files --include-technical 

prefetch SRR8147024 --max-size 420000000000
fasterq-dump SRR8147024 --split-files --include-technical

prefetch SRR8147025 --max-size 420000000000
fasterq-dump SRR8147025 --split-files --include-technical


#concatenate reads
cp SRR8147022_1.fastq Undetermined_S0_L001_R1_001.fastq
gzip Undetermined_S0_L001_R1_001.fastq

cp SRR8147022_2.fastq Undetermined_S0_L001_R2_001.fastq
gzip Undetermined_S0_L001_R2_001.fastq

cp SRR8147022_3.fastq Undetermined_S0_L001_R3_001.fastq
gzip Undetermined_S0_L001_R3_001.fastq

cp SRR8147022_4.fastq Undetermined_S0_L001_R4_001.fastq
gzip Undetermined_S0_L001_R4_001.fastq

cp SRR8147023_1.fastq Undetermined_S0_L002_R1_001.fastq
pigz Undetermined_S0_L002_R1_001.fastq

cp SRR8147023_2.fastq Undetermined_S0_L002_R2_001.fastq
pigz Undetermined_S0_L002_R2_001.fastq

cp SRR8147023_3.fastq Undetermined_S0_L002_R3_001.fastq
pigz Undetermined_S0_L002_R3_001.fastq

cp SRR8147023_4.fastq Undetermined_S0_L002_R4_001.fastq
pigz Undetermined_S0_L002_R4_001.fastq

cp SRR8147024_1.fastq Undetermined_S0_L003_R1_001.fastq
pigz Undetermined_S0_L003_R1_001.fastq

cp SRR8147024_2.fastq Undetermined_S0_L003_R2_001.fastq
pigz Undetermined_S0_L003_R2_001.fastq

cp SRR8147024_3.fastq Undetermined_S0_L003_R3_001.fastq
pigz Undetermined_S0_L003_R3_001.fastq

cp SRR8147024_4.fastq Undetermined_S0_L003_R4_001.fastq
pigz Undetermined_S0_L003_R4_001.fastq

cp SRR8147025_1.fastq Undetermined_S0_L004_R1_001.fastq
pigz Undetermined_S0_L004_R1_001.fastq

cp SRR8147025_2.fastq Undetermined_S0_L004_R2_001.fastq
pigz Undetermined_S0_L004_R2_001.fastq

cp SRR8147025_3.fastq Undetermined_S0_L004_R3_001.fastq
pigz Undetermined_S0_L004_R3_001.fastq

cp SRR8147025_4.fastq Undetermined_S0_L004_R4_001.fastq
pigz Undetermined_S0_L004_R4_001.fastq

