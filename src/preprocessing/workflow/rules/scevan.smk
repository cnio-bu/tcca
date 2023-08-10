rule cna_crc_florian_uhlitz:
    input:
        seurat_list=rules.sc_crc_florian_uhlitz_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/crc_florian_uhlitz_annotated.rds",
    params:
        cna_res=f"{results}/cna/crc_florian_uhlitz"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/crc_florian_uhlitz_scevan.R"

rule cna_oc_ec_matthew_regner:
    input:
        seurat_list=rules.sc_oc_ec_matthew_regner_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/oc_ec_matthew_regner_annotated.rds",
    params:
        cna_res=f"{results}/cna/oc_ec_matthew_regner"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/oc_ec_matthew_regner_scevan.R"

rule cna_nsclc_stefan_salcher:
    input:
        seurat_list=rules.sc_nsclc_stefan_salcher_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/nsclc_stefan_salcher_annotated.rds",
    params:
        cna_res=f"{results}/cna/nsclc_stefan_salcher"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/nsclc_stefan_salcher_scevan.R"

rule cna_esca_xiannian_zhang:
    input:
        seurat_list=rules.sc_esca_xiannian_zhang_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/esca_xiannian_zhang_annotated.rds",
    params:
        cna_res=f"{results}/cna/esca_xiannian_zhang"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/esca_xiannian_zhang_scevan.R"

rule cna_uvm_michael_durante:
    input:
        seurat_list=rules.sc_uvm_michael_durante_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/uvm_michael_durante_annotated.rds",
    params:
        cna_res=f"{results}/cna/uvm_michael_durante"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/uvm_michael_durante_scevan.R"

rule cna_bmets_youmna_kfoury:
    input:
        seurat_list=rules.sc_bmets_youmna_kfoury_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/bmets_youmna_kfoury_annotated.rds",
    params:
        cna_res=f"{results}/cna/bmets_youmna_kfoury"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/bmets_youmna_kfoury_scevan.R"

rule cna_eac_thomas_carroll:
    input:
        seurat_list=rules.sc_eac_thomas_carroll_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/eac_thomas_carroll_annotated.rds",
    params:
        cna_res=f"{results}/cna/eac_thomas_carroll"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/eac_thomas_carroll_scevan.R"

rule cna_cc_xiaosong_lu:
    input:
        seurat_list=rules.sc_cc_xiaosong_lu_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/cc_xiaosong_lu_annotated.rds",
    params:
        cna_res=f"{results}/cna/cc_xiaosong_lu"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/cc_xiaosong_lu_scevan.R"

rule cna_gbm_nourhan_abdelfattah:
    input:
        seurat_list=rules.sc_gbm_nourhan_abdelfattah_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/gbm_nourhan_abdelfattah_annotated.rds",
    params:
        cna_res=f"{results}/cna/gbm_nourhan_abdelfattah"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/gbm_nourhan_abdelfattah_scevan.R"

rule cna_skcm_chao_zhang:
    input:
        seurat_list=rules.sc_skcm_chao_zhang_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/skcm_chao_zhang_annotated.rds",
    params:
        cna_res=f"{results}/cna/skcm_chao_zhang"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/skcm_chao_zhang_scevan.R"

rule cna_brca_bhupinder_pal:
    input:
        seurat_list=rules.sc_brca_bhupinder_pal_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/brca_bhupinder_pal_annotated.rds",
    params:
        cna_res=f"{results}/cna/brca_bhupinder_pal"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/brca_bhupinder_pal_scevan.R"

rule cna_cell_lines_gabriella_kinker:
    input:
        seurat_list=rules.sc_cell_lines_gabriella_kinker_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/cell_lines_gabriella_kinker_annotated.rds",
    params:
        cna_res=f"{results}/cna/cell_lines_gabriella_kinker"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/cell_lines_gabriella_kinker_scevan.R"

rule cna_prad_sujun_chen:
    input:
        seurat_list=rules.sc_prad_sujun_chen_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/prad_sujun_chen_annotated.rds",
    params:
        cna_res=f"{results}/cna/prad_sujun_chen"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/prad_sujun_chen_scevan.R"

rule cna_pdac_junya_peng:
    input:
        seurat_list=rules.sc_pdac_junya_peng_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/pdac_junya_peng_annotated.rds",
    params:
        cna_res=f"{results}/cna/pdac_junya_peng"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/pdac_junya_peng_scevan.R"

rule cna_chol_min_zhang:
    input:
        seurat_list=rules.sc_chol_min_zhang_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/chol_min_zhang_annotated.rds",
    params:
        cna_res=f"{results}/cna/chol_min_zhang"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/chol_min_zhang_scevan.R"

rule cna_mmieloma_stephan_tirier:
    input:
        seurat_list=rules.sc_mmieloma_stephan_tirier_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/mmieloma_stephan_tirier_annotated.rds",
    params:
        cna_res=f"{results}/cna/mmieloma_stephan_tirier"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/mmieloma_stephan_tirier_scevan.R"

rule cna_rcell_kevin_bi:
    input:
        seurat_list=rules.sc_rcell_kevin_bi_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/rcell_kevin_bi_annotated.rds",
    params:
        cna_res=f"{results}/cna/rcell_kevin_bi"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/rcell_kevin_bi_scevan.R"

rule cna_pancancer_junbin_qian:
    input:
        seurat_list=rules.sc_pancancer_junbin_qian_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/pancancer_junbin_qian_annotated.rds",
    params:
        cna_res=f"{results}/cna/pancancer_junbin_qian"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/pancancer_junbin_qian_scevan.R"

rule cna_breast_sunny_wu:
    input:
        seurat_list=rules.sc_breast_sunny_wu_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/breast_sunny_wu_annotated.rds",
    params:
        cna_res=f"{results}/cna/breast_sunny_wu"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/breast_sunny_wu_scevan.R"

rule cna_aml_audrey_lasry:
    input:
        seurat_list=rules.sc_aml_audrey_lasry_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/aml_audrey_lasry_annotated.rds",
    params:
        cna_res=f"{results}/cna/aml_audrey_lasry"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/aml_audrey_lasry_scevan.R"

rule cna_bcc_catherine_dyao:
    input:
        seurat_list=rules.sc_bcc_catherine_dyao_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/bcc_catherine_dyao_annotated.rds",
    params:
        cna_res=f"{results}/cna/bcc_catherine_dyao"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/bcc_catherine_dyao_scevan.R"

rule cna_pancancer_dalia_barkley:
    input:
        seurat_list=rules.sc_pancancer_dalia_barkley_seurat.output.seurat_list,
    output:
        annotated_list=f"{results}/seurat/annotated/pancancer_dalia_barkley_annotated.rds",
    params:
        cna_res=f"{results}/cna/pancancer_dalia_barkley"
    threads: 10
    resources:
        mem_mb=300000,
        walltime=240,
    conda:
        "scevan"
    script:
        "../scripts/pancancer_dalia_barkley_scevan.R"