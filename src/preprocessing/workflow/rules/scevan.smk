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