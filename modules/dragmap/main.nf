nextflow.enable.dsl=2

process Dragmap {

    label 'large_cpu'
    tag 'Preprocessing'

    input: 
        tuple val(meta), path(trimmed_reads)
 
    output:
        tuple val(meta), path("${meta.id}/bam/${meta.id}.sam"), emit: aligned_sam
        // Align report
        path("${meta.id}/bam/${meta.id}.insert-stats.tab"), emit: insert_stat
        path("${meta.id}/bam/${meta.id}.mapping_metrics.csv"), emit: mapping_metrics

    script:
        def rg_args = "--RGID ${meta.flowcell_id} --RGSM ${meta.id}"
        def input_args = trimmed_reads.size() == 1 ? "-1 ${trimmed_reads[0]}" : "-1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]}"

        """
            mkdir -p ${meta.id}/bam
            
            dragen-os \\
                ${rg_args} \\
                --num-threads ${task.cpus} \\
                ${task.ext.args.trim()} \\
                ${input_args} \\
                --output-directory ${meta.id}/bam \\
                --output-file-prefix ${meta.id}
        """
}