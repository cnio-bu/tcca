rule mod_brmets_hugo_gonzalez:
    input:
        bc_list=rules.bc_brmets_hugo_gonzalez.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/brmets_hugo_gonzalez")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_pancancer_dalia_barkley:
    input:
        bc_list=rules.bc_pancancer_dalia_barkley.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/pancancer_dalia_barkley")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_breast_sunny_wu:
    input:
        bc_list=rules.bc_breast_sunny_wu.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/breast_sunny_wu")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_pancancer_sunny_wu:
    input:
        bc_list=rules.bc_pancancer_sunny_wu.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/pancancer_sunny_wu")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_luad_kim_nayoung:
    input:
        bc_list=rules.bc_luad_kim_nayoung.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/luad_kim_nayoung")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_urothelial_chen:
    input:
        bc_list=rules.bc_urothelial_chen.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/urothelial_chen")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_adrenalnb_rui_chong:
    input:
        bc_list=rules.bc_adrenalnb_rui_chong.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/adrenalnb_rui_chong")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_luad_philip_bisschof:
    input:
        bc_list=rules.bc_luad_philip_bisschof.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/luad_philip_bisschof")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_pdac_junya_peng:
    input:
        bc_list=rules.bc_pdac_junya_peng.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/pdac_junya_peng")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_synovial_jerby_arnon:
    input:
        bc_list=rules.bc_synovial_jerby_arnon.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/synovial_jerby_arnon")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_rcell_kevin_bi:
    input:
        bc_list=rules.bc_rcell_kevin_bi.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/rcell_kevin_bi")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_pancancer_junbin_qian:
    input:
        bc_list=rules.bc_pancancer_junbin_qian.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/pancancer_junbin_qian")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_rcell_r_li:
    input:
        bc_list=rules.bc_rcell_r_li.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/rcell_r_li")
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_cll_ramon_massoni:
    input:
        bc_list=rules.bc_cll_ramon_massoni.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/cll_ramon_massoni")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_aml_audrey_lasry:
    input:
        bc_list=rules.bc_aml_audrey_lasry.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/aml_audrey_lasry")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_gbm_nourhan_abdelfattah:
    input:
        bc_list=rules.bc_gbm_nourhan_abdelfattah.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/gbm_nourhan_abdelfattah")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_cc_xiaosong_lu:
    input:
        bc_list=rules.bc_cc_xiaosong_lu.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/cc_xiaosong_lu")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_crc_florian_uhlitz:
    input:
        bc_list=rules.bc_crc_florian_uhlitz.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/crc_florian_uhlitz")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_pleural_rui_dong:
    input:
        bc_list=rules.bc_pleural_rui_dong.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/pleural_rui_dong")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_bone_yun_liu:
    input:
        bc_list=rules.bc_bone_yun_liu.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/bone_yun_liu")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_mmieloma_stephan_tirier:
    input:
        bc_list=rules.bc_mmieloma_stephan_tirier.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/mmieloma_stephan_tirier")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_all_maxime_caron:
    input:
        bc_list=rules.bc_all_maxime_caron.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/all_maxime_caron")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_bcc_catherine_dyao:
    input:
        bc_list=rules.bc_bcc_catherine_dyao.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/bcc_catherine_dyao")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_uvm_michael_durante:
    input:
        bc_list=rules.bc_uvm_michael_durante.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/uvm_michael_durante")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_esca_xiannian_zhang:
    input:
        bc_list=rules.bc_esca_xiannian_zhang.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/esca_xiannian_zhang")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_brca_bhupinder_pal:
    input:
        bc_list=rules.bc_brca_bhupinder_pal.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/brca_bhupinder_pal")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_prad_sujun_chen:
    input:
        bc_list=rules.bc_prad_sujun_chen.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/prad_sujun_chen")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"


rule mod_chol_min_zhang:
    input:
        bc_list=rules.bc_chol_min_zhang.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/chol_min_zhang")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod_skcm_chao_zhang:
    input:
        bc_list=rules.bc_skcm_chao_zhang.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/skcm_chao_zhang")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod_brmets_jana_biermann:
    input:
        bc_list=rules.bc_brmets_jana_biermann.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/brmets_jana_biermann")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180,
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod_cell_lines_gabriella_kinker:
    input:
        bc_list=rules.bc_cell_lines_gabriella_kinker.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/cell_lines_gabriella_kinker")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=180
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod_bmets_youmna_kfoury:
    input:
        bc_list=rules.bc_bmets_youmna_kfoury.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/bmets_youmna_kfoury")
    threads: get_resource("default_bc", "threads")
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

        
rule mod_nsclc_stefan_salcher:
    input:
        bc_list=rules.bc_nsclc_stefan_salcher.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/nsclc_stefan_salcher")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=1440000,
        walltime=180
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"

rule mod_eac_thomas_carroll:
    input:
        bc_list=rules.bc_eac_thomas_carroll.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/eac_thomas_carroll")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=60
    conda:
        "../envs/drug_modules.yaml"
    script:
      "../scripts/therapeutic_module_extraction.R"

rule mod_oc_ec_matthew_regner:
    input:
        bc_list=rules.bc_oc_ec_matthew_regner.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/oc_ec_matthew_regner")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=60
    conda:
      "../envs/drug_modules.yaml"
    script:
      "../scripts/therapeutic_module_extraction.R"

rule mod_pdac_shu_zhang:
    input:
        bc_list=rules.bc_pdac_shu_zhang.output.bc_list,
    output:
        module_dir=directory(f"{results}/modules/pdac_shu_zhang")
    threads: get_resource("default_bc", "threads"),
    resources:
        mem_mb=get_resource("default_bc", "mem_mb"),
        walltime=get_resource("default_bc", "walltime"),
    conda:
        "../envs/drug_modules.yaml"
    script:
        "../scripts/therapeutic_module_extraction.R"
