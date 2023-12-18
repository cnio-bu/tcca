rule sc_brmets_hugo_gonzalez_seurat:
    input:
        mats=get_brmets_mats,
        features=get_brmets_features,
        barcodes=get_brmets_barcodes,
        annotations=get_brmets_annotation,
    output:
        seurat_list=f"{results}/seurat/raw/brmets_hugo_gonzalez.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/brmets_hugo_gonzalez_seurat.R"


rule sc_pancancer_dalia_barkley_seurat:
    input:
        object_list=f"{raw_data}/pancancer_dalia_barkley/srt.list.primary.all.RData",
    output:
        seurat_list=f"{results}/seurat/raw/pancancer_dalia_barkley.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/pancancer_dalia_barkley_seurat.R"


rule sc_breast_sunny_wu_seurat:
    input:
        mat=f"{raw_data}/breast_sunny_wu/matrix.mtx.gz",
        barcodes=f"{raw_data}/breast_sunny_wu/barcodes.tsv.gz",
        features=f"{raw_data}/breast_sunny_wu/features.tsv.gz",
        metadata=f"{raw_data}/breast_sunny_wu/metadata.csv",
    output:
        seurat_list=f"{results}/seurat/raw/breast_sunny_wu.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    params:
        data_dir=f"{raw_data}/breast_sunny_wu",
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/breast_sunny_wu_seurat.R"


rule sc_pancancer_sunny_wu_seurat:
    input:
        mat=f"{raw_data}/pancancer_sunny_wu/matrix.mtx.gz",
        barcodes=f"{raw_data}/pancancer_sunny_wu/barcodes.tsv.gz",
        features=f"{raw_data}/pancancer_sunny_wu/features.tsv.gz",
        metadata=f"{raw_data}/pancancer_sunny_wu/Wu_etal_2021_metadata.txt",
    output:
        seurat_list=f"{results}/seurat/raw/pancancer_sunny_wu.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    params:
        data_dir=f"{raw_data}/pancancer_sunny_wu",
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/pancancer_sunny_wu_seurat.R"


rule sc_luad_kim_nayoung_seurat:
    input:
        mat_object=f"{raw_data}/luad_kim_nayoung/GSE131907_Lung_Cancer_raw_UMI_matrix.rds",
        annotations=f"{raw_data}/luad_kim_nayoung/GSE131907_Lung_Cancer_cell_annotation.txt",
    output:
        seurat_list=f"{results}/seurat/raw/luad_kim_nayoung.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=144000,
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/luad_kim_nayoung_seurat.R"


rule sc_adrenalnb_rui_chong_seurat:
    input:
        object_list=f"{raw_data}/adrenalnb_rui_chong/human_NB_subset_tumor.rda",
        tumor_data=f"{raw_data}/adrenalnb_rui_chong/tumor_dataset_annotation.csv",
        gland_data=f"{raw_data}/adrenalnb_rui_chong/adrenal_gland_annotation.csv",
        reference_gene_annotation="/storage/scratch01/shared/projects/bc-meta/reference/hgnc_gene_with_protein_product_2023-03-22.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/adrenalnb_rui_chong.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/adrenalnb_rui_chong_seurat.R"


rule sc_luad_philip_bisschof_seurat:
    input:
        metadata=f"{raw_data}/luad_philip_bisschof/patients_metadata.xlsx",
        infercnv_scores=f"{raw_data}/luad_philip_bisschof/infercnv_clone_scores_nsclc.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/luad_philip_bisschof.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    params:
        data_dir=f"{raw_data}/luad_philip_bisschof",
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/luad_philip_bisschof_seurat.R"


rule sc_pdac_junya_peng_seurat:
    input:
        object_list=f"{raw_data}/pdac_junya_peng/pdac_junya_peng.RData",
        celltype=f"{raw_data}/pdac_junya_peng/PAAD_CRA001160_CellMetainfo_table.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/pdac_junya_peng.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=144000,
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/pdac_junya_peng_seurat.R"


rule sc_synovial_jerby_arnon_seurat:
    input:
        object_list=f"{raw_data}/synovial_jerby_arnon/seurat_pre-qc.rds",
    output:
        seurat_list=f"{results}/seurat/raw/synovial_jerby_arnon.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/synovial_jerby_arnon_seurat.R"


rule sc_rcell_kevin_bi_seurat:
    input:
        metadata=f"{raw_data}/rcell_kevin_bi/Final_SCP_Metadata.txt",
    output:
        seurat_list=f"{results}/seurat/raw/rcell_kevin_bi.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    params:
        data_dir=f"{raw_data}/rcell_kevin_bi",
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/rcell_kevin_bi_seurat.R"


rule sc_pancancer_junbin_qian_seurat:
    output:
        seurat_list=f"{results}/seurat/raw/pancancer_junbin_qian.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    params:
        data_dir=f"{raw_data}/pancancer_junbin_qian",
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/pancancer_junbin_qian_seurat.R"


rule sc_rcell_r_li_seurat:
    input:
        mat=f"{raw_data}/rcell_r_li/full_mat_annot.tsv",
        annotations=f"{raw_data}/rcell_r_li/annotations.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/rcell_r_li.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=144000,
        walltime=60,
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/rcell_r_li_seurat.R"


rule sc_cll_ramon_massoni_seurat:
    output:
        seurat_list=f"{results}/seurat/raw/cll_ramon_massoni.rds",
    params:
        data_dir=f"{raw_data}/cll_ramon_massoni",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/cll_ramon_massoni_seurat.R"


rule sc_aml_audrey_lasry_seurat:
    input:
        matrix=f"{raw_data}/aml_audrey_lasry/RNA_soupX.mtx",
        features=f"{raw_data}/aml_audrey_lasry/features_RNA_soupX1.csv",
        cells=f"{raw_data}/aml_audrey_lasry/cells_RNA_soupX1.csv",
        metadata=f"{raw_data}/aml_audrey_lasry/metadata_clustering_w_header_upd.csv",
        additional_metadata=f"{raw_data}/aml_audrey_lasry/combined_metadata.csv",
    output:
        seurat_list=f"{results}/seurat/raw/aml_audrey_lasry.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/aml_audrey_lasry_seurat.R"


rule sc_gbm_nourhan_abdelfattah_seurat:
    input:
        metadata=f"{raw_data}/gbm_nourhan_abdelfattah/Meta_Data_GBMatlas.txt",
    output:
        seurat_list=f"{results}/seurat/raw/gbm_nourhan_abdelfattah.rds",
    params:
        data_dir=f"{raw_data}/gbm_nourhan_abdelfattah/GSE182109_RAW",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/gbm_nourhan_abdelfattah_seurat.R"


rule sc_cc_xiaosong_lu_seurat:
    input:
        metadata=f"{raw_data}/cc_xiaosong_lu/metadata_mod.txt",
    output:
        seurat_list=f"{results}/seurat/raw/cc_xiaosong_lu.rds",
    params:
        data_dir=f"{raw_data}/cc_xiaosong_lu",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/cc_xiaosong_lu_seurat.R"


rule sc_pleural_rui_dong_seurat:
    input:
        metadata=f"{raw_data}/pleural_rui_dong/PPB_GSE163678_CellMetainfo_table.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/pleural_rui_dong.rds",
    params:
        data_dir=f"{raw_data}/pleural_rui_dong/GSE163678_RAW"
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/pleural_rui_dong_seurat.R"


rule sc_bone_yun_liu_seurat:
    input:
        metadata=f"{raw_data}/bone_yun_liu/OS_GSE162454_CellMetainfo_table.tsv",
        reference_gene_annotation="/storage/scratch01/shared/projects/bc-meta/reference/hgnc_gene_with_protein_product_2023-03-22.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/bone_yun_liu.rds",
    params:
        data_dir=f"{raw_data}/bone_yun_liu",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/bone_yun_liu_seurat.R"


rule sc_mmieloma_stephan_tirier_seurat:
    input:
        metadata=f"{raw_data}/mmieloma_stephan_tirier/GSE161801_K43R_metadata_table.csv",
    output:
        seurat_list=f"{results}/seurat/raw/mmieloma_stephan_tirier.rds",
    params:
        data_dir=f"{raw_data}/mmieloma_stephan_tirier/GSE161801_RAW",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/mmieloma_stephan_tirier_seurat.R"


rule sc_urothelial_chen_seurat:
    input:
        object_list=f"{raw_data}/urothelial_chen/sc_Chen_BUC.rda",
    output:
        seurat_list=f"{results}/seurat/raw/urothelial_chen.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/urothelial_chen_seurat.R"


rule sc_all_maxime_caron_seurat:
    input:
        metadata=f"{raw_data}/all_maxime_caron/ALL_GSE132509_CellMetainfo_table.tsv",
        reference_gene_annotation="/storage/scratch01/shared/projects/bc-meta/reference/hgnc_gene_with_protein_product_2023-03-22.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/all_maxime_caron.rds",
    params:
        data_dir=f"{raw_data}/all_maxime_caron",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/all_maxime_caron_seurat.R"


rule sc_crc_florian_uhlitz_seurat:
    input:
        metadata=f"{raw_data}/crc_florian_uhlitz/CRC_GSE166555_CellMetainfo_table.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/crc_florian_uhlitz.rds",
    params:
        data_dir=f"{raw_data}/crc_florian_uhlitz/GSE166555_RAW"
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/crc_florian_uhlitz_seurat.R"


rule sc_bcc_catherine_dyao_seurat:
    input:
        metadata=f"{raw_data}/bcc_catherine_dyao/BCC_GSE141526_CellMetainfo_table.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/bcc_catherine_dyao.rds",
    params:
        data_dir=f"{raw_data}/bcc_catherine_dyao/GSE141526_RAW"
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/bcc_catherine_dyao_seurat.R"


rule sc_uvm_michael_durante_seurat:
    input:
        metadata=f"{raw_data}/uvm_michael_durante/UVM_GSE139829_CellMetainfo_table.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/uvm_michael_durante.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    params:
        data_dir = f"{raw_data}/uvm_michael_durante/GSM4147091"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/uvm_michael_durante_seurat.R"


rule sc_esca_xiannian_zhang_seurat:
    input:
        metadata=f"{raw_data}/esca_xiannian_zhang/ESCA_GSE160269_CellMetainfo_table.tsv",
        matrix1=f"{raw_data}/esca_xiannian_zhang/GSE160269/CD45neg_UMIs.txt",
        matrix2=f"{raw_data}/esca_xiannian_zhang/GSE160269/CD45pos_UMIs.txt",
    output:
        seurat_list=f"{results}/seurat/raw/esca_xiannian_zhang.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=84000,
        walltime=120,
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/esca_xiannian_zhang_seurat.R"


rule sc_brca_bhupinder_pal_seurat:
    input:
        metadata=f"{raw_data}/brca_bhupinder_pal/BRCA_GSE161529_CellMetainfo_table.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/brca_bhupinder_pal.rds",
    params:
        data_dir=f"{raw_data}/brca_bhupinder_pal/GSE161529_RAW"
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=120,
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/brca_bhupinder_pal_seurat.R"


rule sc_thyroid_weilin_pu_seurat:
    output:
        seurat_list=f"{results}/seurat/raw/thyroid_weilin_pu.rds",
        out_dir=directory(f"{results}/cna/thyroid_weilin_pu/output"),
    params:
        data_dir=f"{raw_data}/thyroid_weilin_pu/GSE184362_RAW",
        cna_res=f"{results}/cna/thyroid_weilin_pu",
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/thyroid_weilin_pu_seurat.R"


rule sc_prad_sujun_chen_seurat:
    input:
        cell_annot=f"{raw_data}/prad_sujun_chen/PRAD_GSE141445_CellMetainfo_table.tsv",
        raw_matrix=f"{raw_data}/prad_sujun_chen/data.raw.matrix.txt",
    output:
        seurat_list=f"{results}/seurat/raw/prad_sujun_chen.rds",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/prad_sujun_chen_seurat.R"


rule sc_skcm_chao_zhang_seurat:
    input:
        metadata=f"{raw_data}/skcm_chao_zhang/GSE215121_RAW/patient_metadata.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/skcm_chao_zhang.rds",
        out_dir=directory(f"{results}/cna/skcm_chao_zhang/output"),
    params:
        data_dir=f"{raw_data}/skcm_chao_zhang/GSE215121_RAW",
        cna_res=f"{results}/cna/skcm_chao_zhang",
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/skcm_chao_zhang_seurat.R"

rule sc_chol_min_zhang_seurat:
    input:
        metadata=f"{raw_data}/chol_min_zhang/CHOL_GSE138709_CellMetainfo_table.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/chol_min_zhang.rds",
    params:
        data_dir=f"{raw_data}/chol_min_zhang/GSE138709_RAW"
    threads: get_resource("defaults", "threads"),
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/chol_min_zhang_seurat.R"

rule sc_brmets_jana_biermann_seurat:
    input:
        cell_annot=f"{raw_data}/brmets_jana_biermann/GSE200218_sc_sn_metadata.csv",
    output:
        seurat_list=f"{results}/seurat/raw/brmets_jana_biermann.rds",
    params:
        data_dir=f"{raw_data}/brmets_jana_biermann/GSE200218_RAW",
    threads: get_resource("defaults", "threads"),
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/brmets_jana_biermann_seurat.R"

rule sc_cell_lines_gabriella_kinker_seurat:
    input:
        rdata=f"{raw_data}/cell_lines_gabriella_kinker/pancancer_annotated.rds",
    output:
        seurat_list=f"{results}/seurat/raw/cell_lines_gabriella_kinker.rds",
    threads: get_resource("defaults", "threads"),
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/cell_lines_gabriella_kinker_seurat.R"

rule sc_bmets_youmna_kfoury_seurat:
    input:
        metadata=f"{raw_data}/bmets_youmna_kfoury/GSE143791_cell.annotation.human.csv",
    output:
        seurat_list=f"{results}/seurat/raw/bmets_youmna_kfoury.rds",
    params:
        data_dir=f"{raw_data}/bmets_youmna_kfoury/GSE143791_RAW",
    threads: get_resource("defaults", "threads"),
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
        "../envs/seurat.yaml"
    script:
      "../scripts/bmets_youmna_kfoury_seurat.R"

rule sc_nsclc_stefan_salcher_seurat:
    input:
        matrix=f"{raw_data}/nsclc_stefan_salcher/extended_atlas_stefan_salcher.rds",
        reference_gene_annotation="/storage/scratch01/shared/projects/bc-meta/reference/hgnc_gene_with_protein_product_2023-03-22.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/nsclc_stefan_salcher.rds",
    threads: get_resource("defaults", "threads"),
    resources:
        mem_mb=180000,
        walltime=120,
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/nsclc_stefan_salcher_seurat.R"
     
rule sc_eac_thomas_carroll_seurat:
    input:
        matrix=f"{raw_data}/eac_thomas_carroll/sce.pub.Rds",
    output:
        seurat_list=f"{results}/seurat/raw/eac_thomas_carroll.rds",
    threads: get_resource("defaults", "threads"),
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
      "../envs/seurat.yaml"
    script:
      "../scripts/eac_thomas_carroll_seurat.R"

rule sc_oc_ec_matthew_regner_seurat:
    input:
        metadata=f"{raw_data}/oc_ec_matthew_regner/barcode_metadata.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/oc_ec_matthew_regner.rds",
    params:
        data_dir=f"{raw_data}/oc_ec_matthew_regner/GSE173682_RAW",
    threads: get_resource("defaults", "threads"),
    resources:
        mem_mb=get_resource("defaults", "mem_mb"),
        walltime=get_resource("defaults", "walltime"),
    conda:
      "../envs/seurat.yaml"
    script:
      "../scripts/oc_ec_matthew_regner_seurat.R"

rule sc_pdac_shu_zhang_seurat:
    input:
        metadata=f"{raw_data}/pdac_shu_zhang/GSE197177_RAW/patient_metadata.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/pdac_shu_zhang.rds",
        out_dir=directory(f"{results}/cna/pdac_shu_zhang/output"),
    params:
        data_dir=f"{raw_data}/pdac_shu_zhang/GSE197177_RAW/",
        cna_res=f"{results}/cna/pdac_shu_zhang",
    threads: 5
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/pdac_shu_zhang_seurat.R"

rule sc_aml_sander_lambo_seurat:
    input:
        geo_to_samples=f"{raw_data}/aml_sander_lambo/geo_to_sample.tsv",
    output:
        seurat_list=f"{results}/seurat/raw/aml_sander_lambo.rds",
    params:
        data_dir=f"{raw_data}/aml_sander_lambo",
    threads: 1
    resources:
        mem_mb=64000,
        walltime=100,
    conda:
      "../envs/seurat.yaml"
    script:
        "../scripts/aml_sander_lambo_seurat.R"