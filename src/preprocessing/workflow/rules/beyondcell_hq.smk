rule bc2_brmets_hugo_gonzalez:
    input:
        malignant_cells=rules.bc_brmets_hugo_gonzalez.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/brmets_hugo_gonzalez.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_pancancer_dalia_barkley:
    input:
        malignant_cells=rules.bc_pancancer_dalia_barkley.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/pancancer_dalia_barkley.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_breast_sunny_wu:
    input:
        malignant_cells=rules.bc_breast_sunny_wu.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/breast_sunny_wu.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_pancancer_sunny_wu:
    input:
        malignant_cells=rules.bc_pancancer_sunny_wu.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/pancancer_sunny_wu.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_luad_kim_nayoung:
    input:
        malignant_cells=rules.bc_luad_kim_nayoung.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/luad_kim_nayoung.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_urothelial_chen:
    input:
        malignant_cells=rules.bc_urothelial_chen.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/urothelial_chen.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_adrenalnb_rui_chong:
    input:
        malignant_cells=rules.bc_adrenalnb_rui_chong.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/adrenalnb_rui_chong.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_luad_philip_bisschof:
    input:
        malignant_cells=rules.bc_luad_philip_bisschof.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/luad_philip_bisschof.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_pdac_junya_peng:
    input:
        malignant_cells=rules.bc_pdac_junya_peng.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/pdac_junya_peng.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_synovial_jerby_arnon:
    input:
        malignant_cells=rules.bc_synovial_jerby_arnon.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/synovial_jerby_arnon.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_rcell_kevin_bi:
    input:
        malignant_cells=rules.bc_rcell_kevin_bi.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/rcell_kevin_bi.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_pancancer_junbin_qian:
    input:
        malignant_cells=rules.bc_pancancer_junbin_qian.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/pancancer_junbin_qian.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_rcell_r_li:
    input:
        malignant_cells=rules.bc_rcell_r_li.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/rcell_r_li.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_cll_ramon_massoni:
    input:
        malignant_cells=rules.bc_cll_ramon_massoni.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/cll_ramon_massoni.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_aml_audrey_lasry:
    input:
        malignant_cells=rules.bc_aml_audrey_lasry.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/aml_audrey_lasry.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_gbm_nourhan_abdelfattah:
    input:
        malignant_cells=rules.bc_gbm_nourhan_abdelfattah.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/gbm_nourhan_abdelfattah.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_cc_xiaosong_lu:
    input:
        malignant_cells=rules.bc_cc_xiaosong_lu.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/cc_xiaosong_lu.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_pleural_rui_dong:
    input:
        malignant_cells=rules.bc_pleural_rui_dong.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/pleural_rui_dong.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_bone_yun_liu:
    input:
        malignant_cells=rules.bc_bone_yun_liu.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/bone_yun_liu.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_mmieloma_stephan_tirier:
    input:
        malignant_cells=rules.bc_mmieloma_stephan_tirier.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/mmieloma_stephan_tirier.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_all_maxime_caron:
    input:
        malignant_cells=rules.bc_all_maxime_caron.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/all_maxime_caron.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_bcc_catherine_dyao:
    input:
        malignant_cells=rules.bc_bcc_catherine_dyao.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/bcc_catherine_dyao.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_uvm_michael_durante:
    input:
        malignant_cells=rules.bc_uvm_michael_durante.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/uvm_michael_durante.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_esca_xiannian_zhang:
    input:
        malignant_cells=rules.bc_esca_xiannian_zhang.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/esca_xiannian_zhang.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_crc_florian_uhlitz:
    input:
        malignant_cells=rules.bc_crc_florian_uhlitz.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/crc_florian_uhlitz.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_brca_bhupinder_pal:
    input:
        malignant_cells=rules.bc_brca_bhupinder_pal.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/brca_bhupinder_pal.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_prad_sujun_chen:
    input:
        malignant_cells=rules.bc_prad_sujun_chen.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/prad_sujun_chen.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule bc2_chol_min_zhang:
    input:
        malignant_cells=rules.bc_chol_min_zhang.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/chol_min_zhang.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule bc2_skcm_chao_zhang:
    input:
        malignant_cells=rules.bc_skcm_chao_zhang.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/skcm_chao_zhang.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule bc2_brmets_jana_biermann:
    input:
        malignant_cells=rules.bc_brmets_jana_biermann.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/brmets_jana_biermann.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=get_resource("default_deps", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule bc2_cell_lines_gabriella_kinker:
    input:
        malignant_cells=rules.bc_cell_lines_gabriella_kinker.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/cell_lines_gabriella_kinker.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule bc2_bmets_youmna_kfoury:
    input:
        malignant_cells=rules.bc_bmets_youmna_kfoury.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/bmets_youmna_kfoury.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule bc2_nsclc_stefan_salcher:
    input:
        malignant_cells=rules.bc_nsclc_stefan_salcher.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/nsclc_stefan_salcher.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule bc2_eac_thomas_carroll:
    input:
        malignant_cells=rules.bc_eac_thomas_carroll.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/eac_thomas_carroll.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule bc2_oc_ec_matthew_regner:
    input:
        malignant_cells=rules.bc_oc_ec_matthew_regner.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/oc_ec_matthew_regner.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule bc2_pdac_shu_zhang:
    input:
        malignant_cells=rules.bc_pdac_shu_zhang.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/pdac_shu_zhang.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule bc2_aml_sander_lambo:
    input:
        malignant_cells=rules.bc_aml_sander_lambo.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/drug_signatures_fold.gmt",
    output:
        bc_list=f"{results}/beyondcell_hq/aml_sander_lambo.rds",
    threads: get_resource("default_deps", "threads")
    resources:
        mem_mb=get_resource("default_deps", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"