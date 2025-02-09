nextflow.enable.dsl=2

process fgbio_ExtractUmisFromBam {

    label 'large_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(unaligned_bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.unaligned.UMIs.bam"), emit: unaligned_umi_bam
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        def read_structure = meta.read_len == 1 ? "--read-structure=${meta.umi_str}" : "--read-structure=${meta.umi_str} ${meta.umi_str}"
        """
            mkdir -p ${meta.output_dir}/bam

            java \\
            -Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -jar /usr/bin/fgbio.jar \\
            --tmp-dir=${params.tmp_dir} \\
            ExtractUmisFromBam \\
            ${read_structure} \\
            ${args} \\
            --input=${unaligned_bam} \\
            --output=${meta.output_dir}/bam/${meta.id}.unaligned.UMIs.bam
        """
}

process fgbio_GroupReadsByUmi {

    label 'large_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(merged_bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.UMI-groupped.bam"), emit: groupped_umi_bam
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        def strategy = meta.read_len == 1 ? "--strategy=edit" : "--strategy=paired"
        """
            mkdir -p ${meta.output_dir}/bam

            java \\
            -Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -jar /usr/bin/fgbio.jar \\
            --tmp-dir=${params.tmp_dir} \\
            GroupReadsByUmi \\
            ${strategy} \\
            ${args} \\
            --input=${merged_bam} \\
            --output=${meta.output_dir}/bam/${meta.id}.UMI-groupped.bam
        """
}

process fgbio_CallDuplexConsensusReads {

    label 'large_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(groupped_umi_bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.consensus.bam"), emit: consensus_bam
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            java \\
            -Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -jar /usr/bin/fgbio.jar \\
            --tmp-dir=${params.tmp_dir} \\
            CallDuplexConsensusReads \\
            ${args} \\
            --input=${groupped_umi_bam} \\
            --output=${meta.output_dir}/bam/${meta.id}.consensus.bam
        """
}

process fgbio_FilterConsensusReads {

    label 'large_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(consensus_bam)
        val ref
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.consensus.filtered.bam"), emit: filtered_consensus_bam
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            java \\
            -Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -jar /usr/bin/fgbio.jar \\
            --tmp-dir=${params.tmp_dir} \\
            FilterConsensusReads \\
            ${args} \\
            --input=${consensus_bam} \\
            --output=${meta.output_dir}/bam/${meta.id}.consensus.filtered.bam \\
            --ref=${ref.fasta}
        """
}



// Metrics
process fgbio_CollectDuplexSeqMetrics {

    label 'standard'
    // tag "${meta.tag}"
    tag 'Metrics'
    
    input: 
        tuple val(meta), path(groupped_umi_bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/metric/${meta.id}.*"), emit: duplexseqmetric
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/metric

            java \\
            -Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -jar /usr/bin/fgbio.jar \\
            --tmp-dir=${params.tmp_dir} \\
            CollectDuplexSeqMetrics \\
            ${args} \\
            --input=${groupped_umi_bam} \\
            --output=${meta.output_dir}/metric/${meta.id}
        """
}


