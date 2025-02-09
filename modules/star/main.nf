nextflow.enable.dsl=2

process Align {

    label 'high'
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(trimmed_reads)
        val ref
    
    output:
        tuple val(meta), path("${meta.id}/bam/${meta.id}.unaligned.bam")
    
    script:
        def rg_string = [
            meta.flowcell_id ? "-RG ${meta.flowcell_id}" : null,
            meta.id ? "-SM ${meta.id}" : null,
            meta.platform ? "-PL ${meta.platform}" : null,
            meta.library ? "-LB ${meta.library}" : null,
            meta.center ? "-CN ${meta.center}" : null,
        ].findAll { it }.join(" ")
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        def read_command = trimmed_reads[0].toString().endsWith('.gz') ? '--readFilesCommand zcat -' : ""
        """
            STAR \\
            --runMode alignReads \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${ref.genome_dir} \\
            ${read_command} \\

        """
}