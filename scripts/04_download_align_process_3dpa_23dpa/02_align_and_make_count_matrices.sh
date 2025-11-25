singularity exec  indrops.sif python /opt/indrops/indrops.py 3dpa.23dpa.yaml filter

singularity exec indrops.sif python /opt/indrops/indrops.py 3dpa.23dpa.yaml identify_abundant_barcodes

#here we need to direct to indrops from this repo (https://github.com/brianjohnhaas/indrops) to run this custom script

singularity exec indrops.sif python cloned/indrops/extract_barcoded_reads.py 3dpa.23dpa.yaml sort

singularity exec indrops.sif python cloned/indrops/extract_barcoded_reads.py 3dpa.23dpa.yaml  quantify --no-bam > 3dpa_23dpa.barcoded_reads.fastq

cloned/indrops/encode_read_number_in_fastq.pl intact_contra.barcoded_reads.fastq > 3dpa_23dpa.barcoded_reads.adj.fastq

singularity exec  indrops.sif bowtie GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic_w_mito.fna -q -p 48 -a --best --strata --chunkmbs 1000 --sam -m 200 -n 1 -l 15 -e 1003dpa_23dpa.barcoded_reads.adj.fastq  >  GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam


#generate per sample SAM files
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > 3dpa1.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep 3dpa1 > 3dpa1.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > 3dpa2.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep 3dpa2 > 3dpa2.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > 3dpa3.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep 3dpa3 > 3dpa3.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > 23dpa4.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep 23dpa4 > 23dpa4.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam| grep '^@' > 23dpa5.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep 23dpa5 > 23dpa5.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > 23dpa6.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep 23dpa6 > 23dpa6.reform.sam 


#now making counts matrices for each
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 3dpa1.reform.sam > NCBI.3dpa1.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 3dpa2.reform.sam > NCBI.3dpa2.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 3dpa3.reform.sam > NCBI.3dpa3.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 23dpa4.reform.sam > NCBI.23dpa4.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 23dpa5.reform.sam > NCBI.23dpa5.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 23dpa6.reform.sam > NCBI.23dpa6.sc.counts.matrix


python process_rownames.py ../outs/20231128_3dpa_23dpa_inDrops/NCBI.23dpa4.sc.counts.matrix ../outs/20231128_3dpa_23dpa_inDrops/NCBI.23dpa4.sc.gene.counts.matrix
python process_rownames.py ../outs/20231128_3dpa_23dpa_inDrops/NCBI.23dpa5.sc.counts.matrix ../outs/20231128_3dpa_23dpa_inDrops/NCBI.23dpa5.sc.gene.counts.matrix
python process_rownames.py ../outs/20231128_3dpa_23dpa_inDrops/NCBI.23dpa6.sc.counts.matrix ../outs/20231128_3dpa_23dpa_inDrops/NCBI.23dpa6.sc.gene.counts.matrix
python process_rownames.py ../outs/20231128_3dpa_23dpa_inDrops/NCBI.3dpa1.sc.counts.matrix ../outs/20231128_3dpa_23dpa_inDrops/NCBI.3dpa1.sc.gene.counts.matrix
python process_rownames.py ../outs/20231128_3dpa_23dpa_inDrops/NCBI.3dpa2.sc.counts.matrix ../outs/20231128_3dpa_23dpa_inDrops/NCBI.3dpa2.sc.gene.counts.matrix
python process_rownames.py ../outs/20231128_3dpa_23dpa_inDrops/NCBI.3dpa3.sc.counts.matrix ../outs/20231128_3dpa_23dpa_inDrops/NCBI.3dpa3.sc.gene.counts.matrix


#count columns for datamash
#awk '{print NF}' NCBI.23dpa4.sc.counts.matrix | sort -nu | tail -n 1 #23896
#awk '{print NF}' NCBI.23dpa5.sc.counts.matrix | sort -nu | tail -n 1 #21939
#awk '{print NF}' NCBI.23dpa6.sc.counts.matrix | sort -nu | tail -n 1 #7687
#awk '{print NF}' NCBI.3dpa1.sc.counts.matrix | sort -nu | tail -n 1 #22811
#awk '{print NF}' NCBI.3dpa2.sc.counts.matrix | sort -nu | tail -n 1 #19231
#awk '{print NF}' NCBI.3dpa3.sc.counts.matrix | sort -nu | tail -n 1 #22071

#collapse isoform to gene
datamash -g 1 -s -H sum 2-23896 < NCBI.23dpa4.sc.gene.counts.matrix > NCBI.23dpa4.sc.summed.gene.counts.matrix
datamash -g 1 -s -H sum 2-21939 < NCBI.23dpa5.sc.gene.counts.matrix > NCBI.23dpa5.sc.summed.gene.counts.matrix
datamash -g 1 -s -H sum 2-7687 < NCBI.23dpa6.sc.gene.counts.matrix > NCBI.23dpa6.sc.summed.gene.counts.matrix
datamash -g 1 -s -H sum 2-22811 < NCBI.3dpa1.sc.gene.counts.matrix > NCBI.3dpa1.sc.summed.gene.counts.matrix
datamash -g 1 -s -H sum 2-19231 < NCBI.3dpa2.sc.gene.counts.matrix > NCBI.3dpa2.sc.summed.gene.counts.matrix
datamash -g 1 -s -H sum 2-22071 < NCBI.3dpa3.sc.gene.counts.matrix > NCBI.3dpa3.sc.summed.gene.counts.matrix
