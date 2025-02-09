nextflow.enable.dsl=2

process MEM {

    label 'large_cpu'
    tag 'Preprocessing'
    // tag "${meta.tag}"
    
    input: 
        tuple val(meta), path(trimmed_reads)
        val ref
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.aligned.sam"), emit: aligned_sam
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            bwa mem \\
            ${args} \\
            -t ${task.cpus} \\
            ${ref.fasta} \\
            ${trimmed_reads[0]} \\
            ${trimmed_reads[1]} \\
            > ${meta.output_dir}/bam/${meta.id}.aligned.sam
        """
}

