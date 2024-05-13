rule mod2_brmets_hugo_gonzalez:
    input:
        bc_list=rules.bc2_brmets_hugo_gonzalez.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/brmets_hugo_gonzalez")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_pancancer_dalia_barkley:
    input:
        bc_list=rules.bc2_pancancer_dalia_barkley.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/pancancer_dalia_barkley")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_breast_sunny_wu:
    input:
        bc_list=rules.bc2_breast_sunny_wu.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/breast_sunny_wu")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_pancancer_sunny_wu:
    input:
        bc_list=rules.bc2_pancancer_sunny_wu.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/pancancer_sunny_wu")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_luad_kim_nayoung:
    input:
        bc_list=rules.bc2_luad_kim_nayoung.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/luad_kim_nayoung")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_urothelial_chen:
    input:
        bc_list=rules.bc2_urothelial_chen.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/urothelial_chen")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_adrenalnb_rui_chong:
    input:
        bc_list=rules.bc2_adrenalnb_rui_chong.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/adrenalnb_rui_chong")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_luad_philip_bisschof:
    input:
        bc_list=rules.bc2_luad_philip_bisschof.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/luad_philip_bisschof")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_pdac_junya_peng:
    input:
        bc_list=rules.bc2_pdac_junya_peng.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/pdac_junya_peng")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_synovial_jerby_arnon:
    input:
        bc_list=rules.bc2_synovial_jerby_arnon.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/synovial_jerby_arnon")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_rcell_kevin_bi:
    input:
        bc_list=rules.bc2_rcell_kevin_bi.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/rcell_kevin_bi")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_pancancer_junbin_qian:
    input:
        bc_list=rules.bc2_pancancer_junbin_qian.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/pancancer_junbin_qian")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_rcell_r_li:
    input:
        bc_list=rules.bc2_rcell_r_li.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/rcell_r_li")
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_cll_ramon_massoni:
    input:
        bc_list=rules.bc2_cll_ramon_massoni.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/cll_ramon_massoni")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_aml_audrey_lasry:
    input:
        bc_list=rules.bc2_aml_audrey_lasry.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/aml_audrey_lasry")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_gbm_nourhan_abdelfattah:
    input:
        bc_list=rules.bc2_gbm_nourhan_abdelfattah.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/gbm_nourhan_abdelfattah")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_cc_xiaosong_lu:
    input:
        bc_list=rules.bc2_cc_xiaosong_lu.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/cc_xiaosong_lu")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_crc_florian_uhlitz:
    input:
        bc_list=rules.bc2_crc_florian_uhlitz.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/crc_florian_uhlitz")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_pleural_rui_dong:
    input:
        bc_list=rules.bc2_pleural_rui_dong.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/pleural_rui_dong")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_bone_yun_liu:
    input:
        bc_list=rules.bc2_bone_yun_liu.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/bone_yun_liu")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_mmieloma_stephan_tirier:
    input:
        bc_list=rules.bc2_mmieloma_stephan_tirier.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/mmieloma_stephan_tirier")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_all_maxime_caron:
    input:
        bc_list=rules.bc2_all_maxime_caron.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/all_maxime_caron")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_bcc_catherine_dyao:
    input:
        bc_list=rules.bc2_bcc_catherine_dyao.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/bcc_catherine_dyao")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_uvm_michael_durante:
    input:
        bc_list=rules.bc2_uvm_michael_durante.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/uvm_michael_durante")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_esca_xiannian_zhang:
    input:
        bc_list=rules.bc2_esca_xiannian_zhang.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/esca_xiannian_zhang")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_brca_bhupinder_pal:
    input:
        bc_list=rules.bc2_brca_bhupinder_pal.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/brca_bhupinder_pal")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_prad_sujun_chen:
    input:
        bc_list=rules.bc2_prad_sujun_chen.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/prad_sujun_chen")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod2_chol_min_zhang:
    input:
        bc_list=rules.bc2_chol_min_zhang.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/chol_min_zhang")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod2_skcm_chao_zhang:
    input:
        bc_list=rules.bc2_skcm_chao_zhang.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/skcm_chao_zhang")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod2_brmets_jana_biermann:
    input:
        bc_list=rules.bc2_brmets_jana_biermann.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/brmets_jana_biermann")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod2_cell_lines_gabriella_kinker:
    input:
        bc_list=rules.bc2_cell_lines_gabriella_kinker.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/cell_lines_gabriella_kinker")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod2_bmets_youmna_kfoury:
    input:
        bc_list=rules.bc2_bmets_youmna_kfoury.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/bmets_youmna_kfoury")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

        
rule mod2_nsclc_stefan_salcher:
    input:
        bc_list=rules.bc2_nsclc_stefan_salcher.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/nsclc_stefan_salcher")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=1440000,
        walltime=180
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod2_eac_thomas_carroll:
    input:
        bc_list=rules.bc2_eac_thomas_carroll.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/eac_thomas_carroll")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=60
    conda:
        "../envs/drug_modules.yaml"
    script:
      "../scripts/therapeutic_module_extraction.R"

rule mod2_oc_ec_matthew_regner:
    input:
        bc_list=rules.bc2_oc_ec_matthew_regner.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/oc_ec_matthew_regner")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=60
    conda:
      "../envs/drug_modules.yaml"
    script:
      "../scripts/therapeutic_module_extraction.R"

rule mod2_pdac_shu_zhang:
    input:
        bc_list=rules.bc2_pdac_shu_zhang.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/pdac_shu_zhang")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod2_aml_sander_lambo:
    input:
        bc_list=rules.bc2_aml_sander_lambo.output.bc_list,
        moas_table= "../../reference/final_moas - Collapsed.tsv",
    output:
        module_dir=directory(f"{results}/modules_hq/aml_sander_lambo")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"
