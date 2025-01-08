RUNNING AND PREPROCESSING SCEVAN FOR BC META

1. scevan_all_samples: compute SCEVAN for all studies
   scevan_all_samples_2ndrun: compute SCEVAN for biggest studies
2. scevan_metadata: get clonality data per study
3. combine_clonality_metadata: combine all clonality data in one table
4. scevan_cellxgene: generate genes CNV table per study
5. filter_cellxgene_tables: filter genes CNV tables to lvl 3
6. tsv_to_bpcells: transform CNV tables to BPCells matrixes
7. final_cellxgene_v5_cnafill: generate 1-layer Seurat v5 object from BPCells matrixes
8. sketch_rds: generate an sketch from the Seurat v5 object