process samtools_Sort {

    label 'standard'
    tag 'Preprocess'
    // tag "${meta.tag}"
    
    input: 
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.${meta.prefix}.sorted.bam"), emit: sort_bam
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            samtools \\
            sort \\
            -@ ${task.cpus} \\
            ${args} \\
            -o ${meta.output_dir}/bam/${meta.id}.${meta.prefix}.sorted.bam \\
            ${bam}
        """
}

process samtools_Indexing {

    label 'standard'
    tag 'Preprocess'
    // tag "${meta.tag}"
    
    input: 
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${bam}.bai"), emit: index
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            samtools \\
            index \\
            -@ ${task.cpus} \\
            ${args} \\
            ${bam}
        """
}

process samtools_Flagstat {

    label 'standard'
    tag 'Metrics'
    // tag "${meta.tag}"
    
    input: 
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/metric/${meta.id}.flagstat.json"), emit: flagstat
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/metric

            samtools \\
            flagstat \\
            -@ ${task.cpus} \\
            ${args} \\
            ${bam} \\
            > ${meta.output_dir}/metric/${meta.id}.flagstat.json
        """
}