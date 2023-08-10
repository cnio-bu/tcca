rule bc_brmets_hugo_gonzalez:
    input:
        seurat_list=rules.sc_brmets_hugo_gonzalez_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/brmets_hugo_gonzalez.rds",
        bc_list=f"{results}/beyondcell/brmets_hugo_gonzalez.rds",
        report=f"{results}/reports/cells_brmets_hugo_gonzalez.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/brmets_hugo_gonzalez_beyondcell.R"


rule bc_pancancer_dalia_barkley:
    input:
        seurat_list=rules.sc_pancancer_dalia_barkley_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/pancancer_dalia_barkley.rds",
        bc_list=f"{results}/beyondcell/pancancer_dalia_barkley.rds",
        report=f"{results}/reports/cells_pancancer_dalia_barkley.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/pancancer_dalia_barkley_beyondcell.R"


rule bc_breast_sunny_wu:
    input:
        seurat_list=rules.sc_breast_sunny_wu_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/breast_sunny_wu.rds",
        bc_list=f"{results}/beyondcell/breast_sunny_wu.rds",
        report=f"{results}/reports/cells_breast_sunny_wu.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/breast_sunny_wu_beyondcell.R"


rule bc_pancancer_sunny_wu:
    input:
        seurat_list=rules.sc_pancancer_sunny_wu_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/pancancer_sunny_wu.rds",
        bc_list=f"{results}/beyondcell/pancancer_sunny_wu.rds",
        report=f"{results}/reports/cells_pancancer_sunny_wu.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/pancancer_sunny_wu_beyondcell.R"


rule bc_luad_kim_nayoung:
    input:
        seurat_list=rules.sc_luad_kim_nayoung_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/luad_kim_nayoung.rds",
        bc_list=f"{results}/beyondcell/luad_kim_nayoung.rds",
        report=f"{results}/reports/cells_luad_kim_nayoung.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/luad_kim_nayoung_beyondcell.R"


rule bc_urothelial_chen:
    input:
        seurat_list=rules.sc_urothelial_chen_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/urothelial_chen.rds",
        bc_list=f"{results}/beyondcell/urothelial_chen.rds",
        report=f"{results}/reports/cells_urothelial_chen.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/urothelial_chen_beyondcell.R"


rule bc_adrenalnb_rui_chong:
    input:
        seurat_list=rules.sc_adrenalnb_rui_chong_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/adrenalnb_rui_chong.rds",
        bc_list=f"{results}/beyondcell/adrenalnb_rui_chong.rds",
        report=f"{results}/reports/cells_adrenalnb_rui_chong.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/adrenalnb_rui_chong_beyondcell.R"


rule bc_luad_philip_bisschof:
    input:
        seurat_list=rules.sc_luad_philip_bisschof_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/luad_philip_bisschof.rds",
        bc_list=f"{results}/beyondcell/luad_philip_bisschof.rds",
        report=f"{results}/reports/cells_luad_philip_bisschof.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/luad_philip_bisschof_beyondcell.R"


rule bc_pdac_junya_peng:
    input:
        seurat_list=rules.cna_pdac_junya_peng.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/pdac_junya_peng.rds",
        bc_list=f"{results}/beyondcell/pdac_junya_peng.rds",
        report=f"{results}/reports/cells_pdac_junya_peng.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/pdac_junya_peng_beyondcell.R"


rule bc_synovial_jerby_arnon:
    input:
        seurat_list=rules.sc_synovial_jerby_arnon_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/synovial_jerby_arnon.rds",
        bc_list=f"{results}/beyondcell/synovial_jerby_arnon.rds",
        report=f"{results}/reports/cells_synovial_jerby_arnon.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/synovial_jerby_arnon_beyondcell.R"


rule bc_rcell_kevin_bi:
    input:
        seurat_list=rules.cna_rcell_kevin_bi.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/rcell_kevin_bi.rds",
        bc_list=f"{results}/beyondcell/rcell_kevin_bi.rds",
        report=f"{results}/reports/cells_rcell_kevin_bi.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/rcell_kevin_bi_beyondcell.R"


rule bc_pancancer_junbin_qian:
    input:
        seurat_list=rules.sc_pancancer_junbin_qian_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/pancancer_junbin_qian.rds",
        bc_list=f"{results}/beyondcell/pancancer_junbin_qian.rds",
        report=f"{results}/reports/cells_pancancer_junbin_qian.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/pancancer_junbin_qian_beyondcell.R"


rule bc_rcell_r_li:
    input:
        seurat_list=rules.sc_rcell_r_li_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/rcell_r_li.rds",
        bc_list=f"{results}/beyondcell/rcell_r_li.rds",
        report=f"{results}/reports/cells_rcell_r_li.tsv",
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/rcell_r_li_beyondcell.R"


rule bc_cll_ramon_massoni:
    input:
        seurat_list=rules.sc_cll_ramon_massoni_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/cll_ramon_massoni.rds",
        bc_list=f"{results}/beyondcell/cll_ramon_massoni.rds",
        report=f"{results}/reports/cells_cll_ramon_massoni.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/cll_ramon_massoni_beyondcell.R"


rule bc_aml_audrey_lasry:
    input:
        seurat_list=rules.sc_aml_audrey_lasry_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/aml_audrey_lasry.rds",
        bc_list=f"{results}/beyondcell/aml_audrey_lasry.rds",
        report=f"{results}/reports/cells_aml_audrey_lasry.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/aml_audrey_lasry_beyondcell.R"


rule bc_gbm_nourhan_abdelfattah:
    input:
        seurat_list=rules.cna_gbm_nourhan_abdelfattah.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/gbm_nourhan_abdelfattah.rds",
        bc_list=f"{results}/beyondcell/gbm_nourhan_abdelfattah.rds",
        report=f"{results}/reports/cells_gbm_nourhan_abdelfattah.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/gbm_nourhan_abdelfattah_beyondcell.R"


rule bc_cc_xiaosong_lu:
    input:
        seurat_list=rules.cna_cc_xiaosong_lu.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/cc_xiaosong_lu.rds",
        bc_list=f"{results}/beyondcell/cc_xiaosong_lu.rds",
        report=f"{results}/reports/cells_cc_xiaosong_lu.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/cc_xiaosong_lu_beyondcell.R"


rule bc_pleural_rui_dong:
    input:
        seurat_list=rules.sc_pleural_rui_dong_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/pleural_rui_dong.rds",
        bc_list=f"{results}/beyondcell/pleural_rui_dong.rds",
        report=f"{results}/reports/cells_pleural_rui_dong.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/pleural_rui_dong_beyondcell.R"


rule bc_bone_yun_liu:
    input:
        seurat_list=rules.sc_bone_yun_liu_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/bone_yun_liu.rds",
        bc_list=f"{results}/beyondcell/bone_yun_liu.rds",
        report=f"{results}/reports/cells_bone_yun_liu.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/bone_yun_liu_beyondcell.R"


rule bc_mmieloma_stephan_tirier:
    input:
        seurat_list=rules.cna_mmieloma_stephan_tirier.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/mmieloma_stephan_tirier.rds",
        bc_list=f"{results}/beyondcell/mmieloma_stephan_tirier.rds",
        report=f"{results}/reports/cells_mmieloma_stephan_tirier.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/mmieloma_stephan_tirier_beyondcell.R"


rule bc_all_maxime_caron:
    input:
        seurat_list=rules.sc_all_maxime_caron_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/all_maxime_caron.rds",
        bc_list=f"{results}/beyondcell/all_maxime_caron.rds",
        report=f"{results}/reports/cells_all_maxime_caron.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/all_maxime_caron_beyondcell.R"

rule bc_crc_florian_uhlitz:
    input:
        seurat_list=rules.cna_crc_florian_uhlitz.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/crc_florian_uhlitz.rds",
        bc_list=f"{results}/beyondcell/crc_florian_uhlitz.rds",
        report=f"{results}/reports/cells_crc_florian_uhlitz.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/crc_florian_uhlitz_beyondcell.R"


rule bc_bcc_catherine_dyao:
    input:
        seurat_list=rules.sc_bcc_catherine_dyao_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/bcc_catherine_dyao.rds",
        bc_list=f"{results}/beyondcell/bcc_catherine_dyao.rds",
        report=f"{results}/reports/cells_bcc_catherine_dyao.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/bcc_catherine_dyao_beyondcell.R"


rule bc_uvm_michael_durante:
    input:
        seurat_list=rules.cna_uvm_michael_durante.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/uvm_michael_durante.rds",
        bc_list=f"{results}/beyondcell/uvm_michael_durante.rds",
        report=f"{results}/reports/cells_uvm_michael_durante.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/uvm_michael_durante_beyondcell.R"


rule bc_esca_xiannian_zhang:
    input:
        seurat_list=rules.cna_esca_xiannian_zhang.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/esca_xiannian_zhang.rds",
        bc_list=f"{results}/beyondcell/esca_xiannian_zhang.rds",
        report=f"{results}/reports/cells_esca_xiannian_zhang.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/esca_xiannian_zhang_beyondcell.R"


rule bc_brca_bhupinder_pal:
    input:
        seurat_list=rules.cna_brca_bhupinder_pal.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/brca_bhupinder_pal.rds",
        bc_list=f"{results}/beyondcell/brca_bhupinder_pal.rds",
        report=f"{results}/reports/cells_brca_bhupinder_pal.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/brca_bhupinder_pal_beyondcell.R"


rule bc_prad_sujun_chen:
    input:
        seurat_list=rules.cna_prad_sujun_chen.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/prad_sujun_chen.rds",
        bc_list=f"{results}/beyondcell/prad_sujun_chen.rds",
        report=f"{results}/reports/cells_prad_sujun_chen.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/prad_sujun_chen_beyondcell.R"


rule bc_chol_min_zhang:
    input:
        seurat_list=rules.cna_chol_min_zhang.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/chol_min_zhang.rds",
        bc_list=f"{results}/beyondcell/chol_min_zhang.rds",
        report=f"{results}/reports/cells_chol_min_zhang.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/chol_min_zhang_beyondcell.R"

rule bc_skcm_chao_zhang:
    input:
        seurat_list=rules.cna_skcm_chao_zhang.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/skcm_chao_zhang.rds",
        bc_list=f"{results}/beyondcell/skcm_chao_zhang.rds",
        report=f"{results}/reports/cells_skcm_chao_zhang.tsv",
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/skcm_chao_zhang_beyondcell.R"

rule bc_brmets_jana_biermann:
    input:
        seurat_list=rules.sc_brmets_jana_biermann_seurat.output.seurat_list,
    output:
        malignant_list=f"{results}/seurat/malignant/brmets_jana_biermann.rds",
        bc_list=f"{results}/beyondcell/brmets_jana_biermann.rds",
        report=f"{results}/reports/cells_brmets_jana_biermann.tsv",
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/brmets_jana_biermann_beyondcell.R"

rule bc_cell_lines_gabriella_kinker:
    input:
        seurat_list=rules.cna_cell_lines_gabriella_kinker.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/cell_lines_gabriella_kinker.rds",
        bc_list=f"{results}/beyondcell/cell_lines_gabriella_kinker.rds",
        report=f"{results}/reports/cells_cell_lines_gabriella_kinker.tsv",
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/cell_lines_gabriella_kinker_beyondcell.R"

rule bc_bmets_youmna_kfoury:
    input:
        seurat_list=rules.cna_bmets_youmna_kfoury.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/bmets_youmna_kfoury.rds",
        bc_list=f"{results}/beyondcell/bmets_youmna_kfoury.rds",
        report=f"{results}/reports/cells_bmets_youmna_kfoury.tsv",
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/bmets_youmna_kfoury_beyondcell.R"

        
rule bc_nsclc_stefan_salcher:
    input:
        seurat_list=rules.cna_nsclc_stefan_salcher.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/nsclc_stefan_salcher.rds",
        bc_list=f"{results}/beyondcell/nsclc_stefan_salcher.rds",
        report=f"{results}/reports/cells_nsclc_stefan_salcher.tsv",
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=1440000,
        walltime=180
    conda:
        "../envs/beyondcell.yaml"
    script:
        "../scripts/nsclc_stefan_salcher_beyondcell.R"

rule bc_eac_thomas_carroll:
    input:
        seurat_list=rules.cna_eac_thomas_carroll.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/eac_thomas_carroll.rds",
        bc_list=f"{results}/beyondcell/eac_thomas_carroll.rds",
        report=f"{results}/reports/cells_eac_thomas_carroll.tsv",
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=60
    conda:
        "../envs/beyondcell.yaml"
    script:
      "../scripts/eac_thomas_carroll_beyondcell.R"

rule bc_oc_ec_matthew_regner:
    input:
        seurat_list=rules.cna_oc_ec_matthew_regner.output.annotated_list,
    output:
        malignant_list=f"{results}/seurat/malignant/oc_ec_matthew_regner.rds",
        bc_list=f"{results}/beyondcell/oc_ec_matthew_regner.rds",
        report=f"{results}/reports/cells_oc_ec_matthew_regner.tsv",
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=60
    conda:
      "../envs/beyondcell.yaml"
    script:
      "../scripts/oc_ec_matthew_regner_beyondcell.R"
