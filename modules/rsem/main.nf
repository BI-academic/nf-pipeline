
process RSEM_CalculateExpression {

    label 'high'
    tag 'Processing'
    
    input: 
        tuple val(meta), path(transcript_bam), path(rsem_ref)
    
    output:
        tuple (
            val(meta), 
            path("${meta.id}/rsem/${meta.id}.genes.results"),
            path("${meta.id}/rsem/${meta.id}.isoforms.results"),
            path("${meta.id}/rsem/${meta.id}.transcript.bam"),
            path("${meta.id}/rsem/${meta.id}.stat/*")
        )

    script:        
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.id}/rsem

            rsem-calculate-expression \\
            -p ${task.cpus} \\
            ${args} \\
            ${transcript_bam} \\
            ${rsem_ref}/gencode_${meta.gtf_ver} \\
            ${meta.id}/rsem/${meta.id}
        """
}