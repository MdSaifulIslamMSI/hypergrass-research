// Nextflow DSL2 module for HyperGraSS / EIMDL v3 (.hgs2 / .hgs3 + Compressed-Domain Genomics)

process EIMDL_COMPRESS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../bioconda/eimdl/meta.yaml"
    container "ghcr.io/mdsaifulislammsi/eimdl:3.0.0"

    input:
    tuple val(meta), path(target)
    path panel
    val mode_flag

    output:
    tuple val(meta), path("*.hgs2"), emit: archive

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def panel_opt = panel ? "--panel ${panel}" : ""
    def mode_opt = mode_flag ?: (panel ? "" : "--collection")
    """
    eimdl compress \\
        ${mode_opt} \\
        ${panel_opt} \\
        --target ${target} \\
        --out ${prefix}.hgs2 \\
        --threads ${task.cpus}
    """
}

process EIMDL_COMPRESS_OMICS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../bioconda/eimdl/meta.yaml"
    container "ghcr.io/mdsaifulislammsi/eimdl:3.0.0"

    input:
    tuple val(meta), path(reads)
    path panel

    output:
    tuple val(meta), path("*.hgs3"), emit: archive

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def panel_opt = panel ? "--panel ${panel}" : ""
    """
    eimdl compress-omics \\
        ${panel_opt} \\
        --target ${reads} \\
        --cds \\
        --out ${prefix}.hgs3 \\
        --threads ${task.cpus}
    """
}

process EIMDL_ANALYZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../../bioconda/eimdl/meta.yaml"
    container "ghcr.io/mdsaifulislammsi/eimdl:3.0.0"

    input:
    tuple val(meta), path(archive)
    path panel

    output:
    tuple val(meta), path("*.analyze.json"), emit: report

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def panel_opt = panel ? "--panel ${panel}" : ""
    """
    eimdl analyze \\
        --archive ${archive} \\
        ${panel_opt} \\
        --threads ${task.cpus} \\
        --json > ${prefix}.analyze.json
    """
}

process EIMDL_VARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../../bioconda/eimdl/meta.yaml"
    container "ghcr.io/mdsaifulislammsi/eimdl:3.0.0"

    input:
    tuple val(meta), path(archive)
    path panel

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.variants.json"), emit: summary

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    eimdl variants \\
        --archive ${archive} \\
        --panel ${panel} \\
        --out ${prefix}.vcf \\
        --json > ${prefix}.variants.json
    """
}

process EIMDL_KMER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../../bioconda/eimdl/meta.yaml"
    container "ghcr.io/mdsaifulislammsi/eimdl:3.0.0"

    input:
    tuple val(meta), path(archive)
    path panel
    val query_seq

    output:
    tuple val(meta), path("*.kmer.json"), emit: hits

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def panel_opt = panel ? "--panel ${panel}" : ""
    """
    eimdl kmer \\
        --archive ${archive} \\
        ${panel_opt} \\
        --query ${query_seq} \\
        --both-strands \\
        --json > ${prefix}.kmer.json
    """
}
