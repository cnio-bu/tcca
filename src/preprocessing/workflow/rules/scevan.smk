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