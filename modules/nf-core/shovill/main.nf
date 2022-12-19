process SHOVILL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::shovill=1.1.0" : null) 
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shovill%3A1.1.0--hdfd78af_1' :
        'quay.io/biocontainers/shovill:1.1.0--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("$meta.id/*.contigs.fa")                       , emit: contigs
    tuple val(meta), path("$meta.id/shovill.corrections")                , emit: corrections
    tuple val(meta), path("$meta.id/shovill.log")                        , emit: log
    tuple val(meta), path("$meta.id/{skesa,spades,megahit,velvet}.fasta"), emit: raw_contigs
    tuple val(meta), path("$meta.id/contigs.{fastg,gfa,LastGraph}")      , optional:true, emit: gfa
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def memory = task.memory.toGiga()
    """
    mkdir -p $meta.id
    shovill \\
        --R1 ${reads[0]} \\
        --R2 ${reads[1]} \\
        $args \\
        --cpus $task.cpus \\
        --ram $memory \\
        --outdir ./$meta.id \\
        --force

    for $file in ${meta.id}/*; do
    mv ${file} ${meta.id}.${file}
    done
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: \$(echo \$(shovill --version 2>&1) | sed 's/^.*shovill //')
    END_VERSIONS
    """
}
