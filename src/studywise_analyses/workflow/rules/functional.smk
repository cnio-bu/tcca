rule fc_brmets_hugo_gonzalez:
    input:
        malignant_cells=rules.bc_brmets_hugo_gonzalez.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/brmets_hugo_gonzalez.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_pancancer_dalia_barkley:
    input:
        malignant_cells=rules.bc_pancancer_dalia_barkley.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/pancancer_dalia_barkley.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_breast_sunny_wu:
    input:
        malignant_cells=rules.bc_breast_sunny_wu.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/breast_sunny_wu.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_pancancer_sunny_wu:
    input:
        malignant_cells=rules.bc_pancancer_sunny_wu.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/pancancer_sunny_wu.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_luad_kim_nayoung:
    input:
        malignant_cells=rules.bc_luad_kim_nayoung.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/luad_kim_nayoung.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_urothelial_chen:
    input:
        malignant_cells=rules.bc_urothelial_chen.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/urothelial_chen.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_adrenalnb_rui_chong:
    input:
        malignant_cells=rules.bc_adrenalnb_rui_chong.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/adrenalnb_rui_chong.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_luad_philip_bisschof:
    input:
        malignant_cells=rules.bc_luad_philip_bisschof.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/luad_philip_bisschof.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_pdac_junya_peng:
    input:
        malignant_cells=rules.bc_pdac_junya_peng.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/pdac_junya_peng.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_synovial_jerby_arnon:
    input:
        malignant_cells=rules.bc_synovial_jerby_arnon.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/synovial_jerby_arnon.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_rcell_kevin_bi:
    input:
        malignant_cells=rules.bc_rcell_kevin_bi.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/rcell_kevin_bi.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_pancancer_junbin_qian:
    input:
        malignant_cells=rules.bc_pancancer_junbin_qian.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/pancancer_junbin_qian.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_rcell_r_li:
    input:
        malignant_cells=rules.bc_rcell_r_li.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/rcell_r_li.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_cll_ramon_massoni:
    input:
        malignant_cells=rules.bc_cll_ramon_massoni.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/cll_ramon_massoni.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_aml_audrey_lasry:
    input:
        malignant_cells=rules.bc_aml_audrey_lasry.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/aml_audrey_lasry.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_gbm_nourhan_abdelfattah:
    input:
        malignant_cells=rules.bc_gbm_nourhan_abdelfattah.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/gbm_nourhan_abdelfattah.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_cc_xiaosong_lu:
    input:
        malignant_cells=rules.bc_cc_xiaosong_lu.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/cc_xiaosong_lu.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_pleural_rui_dong:
    input:
        malignant_cells=rules.bc_pleural_rui_dong.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/pleural_rui_dong.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_bone_yun_liu:
    input:
        malignant_cells=rules.bc_bone_yun_liu.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/bone_yun_liu.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_mmieloma_stephan_tirier:
    input:
        malignant_cells=rules.bc_mmieloma_stephan_tirier.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/mmieloma_stephan_tirier.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_all_maxime_caron:
    input:
        malignant_cells=rules.bc_all_maxime_caron.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/all_maxime_caron.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_florian_uhlitz:
    input:
        malignant_cells=rules.bc_crc_florian_uhlitz.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/florian_uhlitz.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_bcc_catherine_dyao:
    input:
        malignant_cells=rules.bc_bcc_catherine_dyao.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/bcc_catherine_dyao.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_uvm_michael_durante:
    input:
        malignant_cells=rules.bc_uvm_michael_durante.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/uvm_michael_durante.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_esca_xiannian_zhang:
    input:
        malignant_cells=rules.bc_esca_xiannian_zhang.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/esca_xiannian_zhang.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_crc_florian_uhlitz:
    input:
        malignant_cells=rules.bc_crc_florian_uhlitz.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/crc_florian_uhlitz.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_brca_bhupinder_pal:
    input:
        malignant_cells=rules.bc_brca_bhupinder_pal.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/brca_bhupinder_pal.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_prad_sujun_chen:
    input:
        malignant_cells=rules.bc_prad_sujun_chen.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/prad_sujun_chen.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"


rule fc_chol_min_zhang:
    input:
        malignant_cells=rules.bc_chol_min_zhang.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/chol_min_zhang.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_skcm_chao_zhang:
    input:
        malignant_cells=rules.bc_skcm_chao_zhang.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/skcm_chao_zhang.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_brmets_jana_biermann:
    input:
        malignant_cells=rules.bc_brmets_jana_biermann.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/brmets_jana_biermann.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_cell_lines_gabriella_kinker:
    input:
        malignant_cells=rules.bc_cell_lines_gabriella_kinker.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/cell_lines_gabriella_kinker.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_bmets_youmna_kfoury:
    input:
        malignant_cells=rules.bc_bmets_youmna_kfoury.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/bmets_youmna_kfoury.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_nsclc_stefan_salcher:
    input:
        malignant_cells=rules.bc_nsclc_stefan_salcher.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/nsclc_stefan_salcher.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_eac_thomas_carroll:
    input:
        malignant_cells=rules.bc_eac_thomas_carroll.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/eac_thomas_carroll.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_oc_ec_matthew_regner:
    input:
        malignant_cells=rules.bc_oc_ec_matthew_regner.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/oc_ec_matthew_regner.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_pdac_shu_zhang:
    input:
        malignant_cells=rules.bc_pdac_shu_zhang.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/pdac_shu_zhang.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"

rule fc_aml_sander_lambo:
    input:
        malignant_cells=rules.bc_aml_sander_lambo.output.malignant_list,
        gsets="/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt",
    output:
        bc_list=f"{results}/functional/aml_sander_lambo.rds",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/general_functional_enrichment.R"