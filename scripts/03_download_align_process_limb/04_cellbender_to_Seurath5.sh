#module load Anaconda3/2023.07-2

ptrepack --complevel 5 Intact1.reform_filtered.h5:/matrix h5Seurat_reform/Intact1_so_output_filtered_seurat.h5:/matrix
ptrepack --complevel 5 Intact2.reform_filtered.h5:/matrix h5Seurat_reform/Intact2_so_output_filtered_seurat.h5:/matrix
ptrepack --complevel 5 Contra1.reform_filtered.h5:/matrix h5Seurat_reform/Contra1_so_output_filtered_seurat.h5:/matrix
ptrepack --complevel 5 Contra2.reform_filtered.h5:/matrix h5Seurat_reform/Contra2_so_output_filtered_seurat.h5:/matrix
ptrepack --complevel 5 Contra3.reform_filtered.h5:/matrix h5Seurat_reform/Contra3_so_output_filtered_seurat.h5:/matrix

ptrepack --complevel 5 Intact1.summed_filtered.h5:/matrix h5Seurat_summed/Intact1_so_output_filtered_seurat.h5:/matrix
ptrepack --complevel 5 Intact2.summed_filtered.h5:/matrix h5Seurat_summed/Intact2_so_output_filtered_seurat.h5:/matrix
ptrepack --complevel 5 Contra1.summed_filtered.h5:/matrix h5Seurat_summed/Contra1_so_output_filtered_seurat.h5:/matrix
ptrepack --complevel 5 Contra2.summed_filtered.h5:/matrix h5Seurat_summed/Contra2_so_output_filtered_seurat.h5:/matrix
ptrepack --complevel 5 Contra3.summed_filtered.h5:/matrix h5Seurat_summed/Contra3_so_output_filtered_seurat.h5:/matrix
