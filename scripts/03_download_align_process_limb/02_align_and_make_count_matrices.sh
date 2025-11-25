singularity exec  indrops.sif python /opt/indrops/indrops.py limb.yaml filter

singularity exec indrops.sif python /opt/indrops/indrops.py limb.yaml identify_abundant_barcodes

singularity exec indrops.sif python cloned/indrops/extract_barcoded_reads.py limb.yaml sort

#here we need to direct to indrops from this repo (https://github.com/brianjohnhaas/indrops) to run this custom script
singularity exec indrops.sif python cloned/indrops/extract_barcoded_reads.py Intact_contra.yaml  quantify --no-bam > intact_contra.barcoded_reads.fastq

cloned/indrops/encode_read_number_in_fastq.pl intact_contra.barcoded_reads.fastq > intact_contra.barcoded_reads.adj.fastq

singularity exec --bind indrops.sif bowtie GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.fna  -q -p 48 -a --best --strata --chunkmbs 1000 --sam -m 200 -n 1 -l 15 -e 100 intact_contra.barcoded_reads.adj.fastq  >  GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam

#using  SAMtools/1.16.1 pull out sam files for each sample

samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > Intact1.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep Intact1 > Intact1.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > Intact2.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep Intact1 > Intact1.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > Contra1.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep Contra1 > Contra1.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > Contra2.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep Contra2 > Contra2.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > Contra3.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep Contra3 > Contra3.reform.sam 


#now making counts matrices for each
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam Intact1.reform.sam > NCBI.Intact1.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam Intact2.reform.sam > NCBI.Intact2.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam Contra1.reform.sam > NCBI.Contra1.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam Contra2.reform.sam > NCBI.Contra2.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam Contra3.reform.sam > NCBI.Contra3.sc.counts.matrix


#Intact 1
python process_rownames.py NCBI.Intact1.sc.counts.matrix NCBI.Intact1.sc.gene.counts.matrix
#check number of columns
awk '{print NF}' NCBI.Intact1.sc.gene.counts.matrix | sort -nu | tail -n 1 
#10817 
#combine all isoforms into gene level counts
datamash -g 1 -s -H sum 2-10817 < NCBI.Intact1.sc.gene.counts.matrix > NCBI.Intact1.sc.summed.gene.counts.matrix

#Intact 2
python process_rownames.py NCBI.Intact2.sc.counts.matrix NCBI.Intact2.sc.gene.counts.matrix
awk '{print NF}' NCBI.Intact2.sc.gene.counts.matrix | sort -nu | tail -n 1 
#14099
#combine all isoforms into gene level counts
datamash -g 1 -s -H sum 2-14099 < NCBI.Intact2.sc.gene.counts.matrix > NCBI.Intact2.sc.summed.gene.counts.matrix

#Contra 1
python process_rownames.py NCBI.Contra1.sc.counts.matrix NCBI.Contra1.sc.gene.counts.matrix
#awk '{print NF}' NCBI.Contra1.sc.counts.matrix | sort -nu | tail -n 1 
datamash -g 1 -s -H sum 2-19905 < NCBI.Contra1.sc.gene.counts.matrix > NCBI.Contra1.sc.summed.gene.counts.matrix

#Contra 2
python process_rownames.py NCBI.Contra2.sc.counts.matrix NCBI.Contra2.sc.gene.counts.matrix
#awk '{print NF}' NCBI.Contra2.sc.counts.matrix | sort -nu | tail -n 1 
#4400
datamash -g 1 -s -H sum 2-4400 < NCBI.Contra2.sc.gene.counts.matrix > NCBI.Contra2.sc.summed.gene.counts.matrix

#Contra 3
python process_rownames.py ../outs/20231128_Intact_Contra/NCBI.Contra3.sc.counts.matrix ../outs/20231128_Intact_Contra/NCBI.Contra3.sc.gene.counts.matrix
#awk '{print NF}' NCBI.Contra3.sc.counts.matrix | sort -nu | tail -n 1 
#4870
datamash -g 1 -s -H sum 2-4870 < NCBI.Contra3.sc.gene.counts.matrix > NCBI.Contra3.sc.summed.gene.counts.matrix

