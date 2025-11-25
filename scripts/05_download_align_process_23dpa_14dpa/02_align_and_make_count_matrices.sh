singularity exec  indrops.sif python /opt/indrops/indrops.py 23dpa_14dpa.yaml filter

singularity exec indrops.sif python /opt/indrops/indrops.py 23dpa_14dpa.yaml identify_abundant_barcodes

#here we need to direct to indrops from this repo (https://github.com/brianjohnhaas/indrops) to run this custom script

singularity exec indrops.sif python cloned/indrops/extract_barcoded_reads.py 23dpa_14dpa.yaml sort

singularity exec indrops.sif python cloned/indrops/extract_barcoded_reads.py 23dpa_14dpa.yaml  quantify --no-bam > 23dpa_14dpa.barcoded_reads.fastq

cloned/indrops/encode_read_number_in_fastq.pl intact_contra.barcoded_reads.fastq > 23dpa_14dpa.barcoded_reads.adj.fastq

singularity exec --bind indrops.sif bowtie GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic_w_mito.fna  -q -p 48 -a --best --strata --chunkmbs 1000 --sam -m 200 -n 1 -l 15 -e 100 23dpa_14dpa.barcoded_reads.adj.fastq  >  GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam


#generate per sample SAM files
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > 23dpa1.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep S_1 > 23dpa1.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > 23dpa2.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep S_2 > 23dpa2.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > 23dpa3.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep S_4 > 23dpa3.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep '^@' > 14dpa1.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep S_3 > 14dpa1.reform.sam 
samtools view -h GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam| grep '^@' > 14dpa2.reform.sam
samtools view GCF_040938575.1_UKY_AmexF1_1_rna_from_genomic.reformatted.target.bowtie.adj.sam | grep S_5 > 14dpa2.reform.sam 

#now making counts matrices for each
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 23dpa1.reform.sam > NCBI.23dpa1.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 23dpa2.reform.sam > NCBI.23dpa2.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 23dpa3.reform.sam > NCBI.23dpa3.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 14dpa1.reform.sam > NCBI.14dpa1.sc.counts.matrix
${ScriptDIR}bam_to_count_matrix.pl --max_top_cells 1800000 --min_reads 0 --bam 14dpa2.reform.sam > NCBI.14dpa2.sc.counts.matrix

#cut off NCBI accessions
python process_rownames.py ../outs/20231128_23dpa_14dpa_indrops/NCBI.23dpa1.sc.counts.matrix ../outs/20231128_23dpa_14dpa_indrops/NCBI.23dpa1.sc.gene.counts.matrix
python process_rownames.py ../outs/20231128_23dpa_14dpa_indrops/NCBI.23dpa2.sc.counts.matrix ../outs/20231128_23dpa_14dpa_indrops/NCBI.23dpa2.sc.gene.counts.matrix
python process_rownames.py ../outs/20231128_23dpa_14dpa_indrops/NCBI.23dpa3.sc.counts.matrix ../outs/20231128_23dpa_14dpa_indrops/NCBI.23dpa3.sc.gene.counts.matrix
python process_rownames.py ../outs/20231128_23dpa_14dpa_indrops/NCBI.14dpa1.sc.counts.matrix ../outs/20231128_23dpa_14dpa_indrops/NCBI.14dpa1.sc.gene.counts.matrix
python process_rownames.py ../outs/20231128_23dpa_14dpa_indrops/NCBI.14dpa2.sc.counts.matrix ../outs/20231128_23dpa_14dpa_indrops/NCBI.14dpa2.sc.gene.counts.matrix

cd ../outs/20231128_23dpa_14dpa_indrops/

#awk '{print NF}' NCBI.23dpa1.sc.counts.matrix | sort -nu | tail -n 1 #23308
#awk '{print NF}' NCBI.23dpa2.sc.counts.matrix | sort -nu | tail -n 1 #21937
#awk '{print NF}' NCBI.23dpa3.sc.counts.matrix | sort -nu | tail -n 1 #24705
#awk '{print NF}' NCBI.14dpa1.sc.counts.matrix | sort -nu | tail -n 1 #29319
#awk '{print NF}' NCBI.14dpa2.sc.counts.matrix | sort -nu | tail -n 1 #31735

datamash -g 1 -s -H sum 2-23308 < NCBI.23dpa1.sc.gene.counts.matrix > NCBI.23dpa1.sc.summed.gene.counts.matrix
datamash -g 1 -s -H sum 2-21937 < NCBI.23dpa2.sc.gene.counts.matrix > NCBI.23dpa2.sc.summed.gene.counts.matrix
datamash -g 1 -s -H sum 2-24705 < NCBI.23dpa3.sc.gene.counts.matrix > NCBI.23dpa3.sc.summed.gene.counts.matrix
datamash -g 1 -s -H sum 2-29319 < NCBI.14dpa1.sc.gene.counts.matrix > NCBI.14dpa1.sc.summed.gene.counts.matrix
datamash -g 1 -s -H sum 2-31735 < NCBI.14dpa2.sc.gene.counts.matrix > NCBI.14dpa2.sc.summed.gene.counts.matrix
