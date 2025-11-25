# adaptive_immune_limb_regen

We reprocessed data from [Leigh et al. 2018 ](https://www.nature.com/articles/s41467-018-07604-0). Here you will find all the information on how to reprocess the data and the code we used to re-analyze it. 

We took advantage of this [image](https://hub.docker.com/r/shengqh/indrops)

Which can be obtained via: 

`apptainer pull indrops.sif docker://shengqh/indrops`

# General processing pipeline

## 1. Download, align, and clean
The sample sheet for this data can be above called sample.table.txt. Each set of samples are processed in batches based on how they were sequenced (see sample sheet for clarification). For clarity, we have documented how we downloaded and preprocessed each batch in the scripts directory. These can be found under the directories
01_download_process_align_limb
02_download_process_align_early_and_med_bud
03_download_process_align_WH_and_med_bud

Within each subdirectory you will find the steps we took for download and processing of each sample. We did the following for each. 

1. Download from SRA using SRAtoolkit v3.0.3 via pre-fetch and fasterq-dump
2. Concatenate reads and zip using pigz
3. Using the indrops pipeline described [here](https://github.com/brianjohnhaas/indrops), filter, identify_abundant_barcodes, sort, encode_read_number_in_fastq.pl, and finally align using bowtie. The above noted indrops container in combination with the provided .yaml file in each subdirectory is critical for running this step. 
4. Use SAMtools to generate sample specific SAM files.
5. Generate sample specific count matrices using bam_to_count_matrix.pl (again from [here](https://github.com/brianjohnhaas/indrops)). We intentionally output a very high number of max_top_cells to allow for advanced filtering in Cellbender.
6. Collapse from isoform to genes using process_rownames.py to strip away NBCI accession number and leave gene names and then collapse matrices to gene-level using datamash. 
7. Create h5ad files for each sample in preperation for [Cellbender](https://github.com/broadinstitute/CellBender). Each batch folder contains an Rscript .R file to do so.
8. Run Cellbender
9. Outupt h5ad from Cellbender are then converted to .h5 for use in Seurat.
10. Now these matrices are ready to load into Seurat and following the analysis provided in the code directory. 

## 2. Analyze
