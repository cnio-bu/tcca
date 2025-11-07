rule pc_brmets_hugo_gonzalez:
    input:
        malignant_list=rules.bc_brmets_hugo_gonzalez.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/brmets_hugo_gonzalez")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_pancancer_dalia_barkley:
    input:
        malignant_list=rules.bc_pancancer_dalia_barkley.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/pancancer_dalia_barkley")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_breast_sunny_wu:
    input:
        malignant_list=rules.bc_breast_sunny_wu.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/breast_sunny_wu")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_pancancer_sunny_wu:
    input:
        malignant_list=rules.bc_pancancer_sunny_wu.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/pancancer_sunny_wu")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_luad_kim_nayoung:
    input:
        malignant_list=rules.bc_luad_kim_nayoung.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/luad_kim_nayoung")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_urothelial_chen:
    input:
        malignant_list=rules.bc_urothelial_chen.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/urothelial_chen")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_adrenalnb_rui_chong:
    input:
        malignant_list=rules.bc_adrenalnb_rui_chong.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/adrenalnb_rui_chong")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_luad_philip_bisschof:
    input:
        malignant_list=rules.bc_luad_philip_bisschof.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/luad_philip_bisschof")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_pdac_junya_peng:
    input:
        malignant_list=rules.bc_pdac_junya_peng.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/pdac_junya_peng")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_synovial_jerby_arnon:
    input:
        malignant_list=rules.bc_synovial_jerby_arnon.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/synovial_jerby_arnon")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_rcell_kevin_bi:
    input:
        malignant_list=rules.bc_rcell_kevin_bi.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/rcell_kevin_bi")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_pancancer_junbin_qian:
    input:
        malignant_list=rules.bc_pancancer_junbin_qian.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/pancancer_junbin_qian")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_rcell_r_li:
    input:
        malignant_list=rules.bc_rcell_r_li.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/rcell_r_li")
    threads: get_resource("defaults", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_cll_ramon_massoni:
    input:
        malignant_list=rules.bc_cll_ramon_massoni.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/cll_ramon_massoni")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_aml_audrey_lasry:
    input:
        malignant_list=rules.bc_aml_audrey_lasry.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/aml_audrey_lasry")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_gbm_nourhan_abdelfattah:
    input:
        malignant_list=rules.bc_gbm_nourhan_abdelfattah.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/gbm_nourhan_abdelfattah")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=64000,
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_cc_xiaosong_lu:
    input:
        malignant_list=rules.bc_cc_xiaosong_lu.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/cc_xiaosong_lu")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_pleural_rui_dong:
    input:
        malignant_list=rules.bc_pleural_rui_dong.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/pleural_rui_dong")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_bone_yun_liu:
    input:
        malignant_list=rules.bc_bone_yun_liu.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/bone_yun_liu")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_mmieloma_stephan_tirier:
    input:
        malignant_list=rules.bc_mmieloma_stephan_tirier.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/mmieloma_stephan_tirier")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_all_maxime_caron:
    input:
        malignant_list=rules.bc_all_maxime_caron.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/all_maxime_caron")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"

rule pc_crc_florian_uhlitz:
    input:
        malignant_list=rules.bc_crc_florian_uhlitz.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/crc_florian_uhlitz")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_bcc_catherine_dyao:
    input:
        malignant_list=rules.bc_bcc_catherine_dyao.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/bcc_catherine_dyao")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_uvm_michael_durante:
    input:
        malignant_list=rules.bc_uvm_michael_durante.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/uvm_michael_durante")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_esca_xiannian_zhang:
    input:
        malignant_list=rules.bc_esca_xiannian_zhang.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/esca_xiannian_zhang")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=180,
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_brca_bhupinder_pal:
    input:
        malignant_list=rules.bc_brca_bhupinder_pal.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/brca_bhupinder_pal")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=180,
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_prad_sujun_chen:
    input:
        malignant_list=rules.bc_prad_sujun_chen.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/prad_sujun_chen")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_chol_min_zhang:
    input:
        malignant_list=rules.bc_chol_min_zhang.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/chol_min_zhang")
    threads: get_resource("default_pc", "threads")
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"


rule pc_skcm_chao_zhang:
    input:
        malignant_list=rules.bc_skcm_chao_zhang.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/skcm_chao_zhang")
    threads: get_resource("default_pc", "threads"),
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"

rule pc_brmets_jana_biermann:
    input:
        malignant_list=rules.bc_brmets_jana_biermann.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/brmets_jana_biermann")
    threads: get_resource("default_pc", "threads"),
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"

rule pc_cell_lines_gabriella_kinker:
    input:
        malignant_list=rules.bc_cell_lines_gabriella_kinker.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/cell_lines_gabriella_kinker")
    threads: get_resource("default_pc", "threads"),
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=180
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"

       
rule pc_nsclc_stefan_salcher:
    input:
        malignant_list=rules.bc_nsclc_stefan_salcher.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/nsclc_stefan_salcher")
    threads: get_resource("default_pc", "threads"),
    resources:
        mem_mb=1440000,
        runtime=180
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"

rule pc_eac_thomas_carroll:
    input:
        malignant_list=rules.bc_eac_thomas_carroll.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/eac_thomas_carroll")
    threads: get_resource("default_pc", "threads"),
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=60
    conda:
        "../envs/perception.yaml"
    script:
      "../scripts/general_perception.R"

rule pc_oc_ec_matthew_regner:
    input:
        malignant_list=rules.bc_oc_ec_matthew_regner.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/oc_ec_matthew_regner")
    threads: get_resource("default_pc", "threads"),
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=60
    conda:
      "../envs/perception.yaml"
    script:
      "../scripts/general_perception.R"

rule pc_pdac_shu_zhang:
    input:
        malignant_list=rules.bc_pdac_shu_zhang.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/pdac_shu_zhang")
    threads: get_resource("default_pc", "threads"),
    resources:
        mem_mb=get_resource("default_pc", "mem_mb"),
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"

rule pc_aml_sander_lambo:
    input:
        malignant_list=rules.bc_aml_sander_lambo.output.malignant_list,
        drug_models=f"{results}/perception/drug_models/FDA_approved_drugs_models.rds",
    output: perception_mat=directory(f"{results}/perception/v5/aml_sander_lambo")
    threads: get_resource("default_pc", "threads"),
    resources:
        mem_mb=180000,
        runtime=get_resource("default_pc", "runtime"),
    conda:
        "../envs/perception.yaml"
    script:
        "../scripts/general_perception.R"