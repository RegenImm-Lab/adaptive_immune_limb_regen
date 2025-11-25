apptainer build Seurat.sif docker://satijalab/seurat:v4.3.0

apptainer  exec  ${ScriptDIR}Seurat.sif Rscript Rscript.contra.intact.save.R
