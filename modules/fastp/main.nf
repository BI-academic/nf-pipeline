nextflow.enable.dsl=2

process Fastp {

    label 'standard'
    tag 'Preprocessing'
    // tag "${meta.tag}"

    input: 
        tuple val(meta), path(reads)
 
    output:
        tuple val(meta), path("${meta.output_dir}/fastq/*.trimmed_r{1,2}.fastq.gz"), emit: trimmed_reads
        // Trimming report
        path "${meta.output_dir}/fastq/*.json", emit: json
        path "${meta.output_dir}/fastq/*.html", emit: html
        // path '${sample_id}/fastp/*.log', emit: log
        path "${meta.output_dir}/fastq/*.failed.fastq.gz", optional: true, emit: failed_out
        path "${meta.output_dir}/fastq/*.unpaired_*.fastq.gz", optional: true, emit: unpaired_reads

    script:
        def input_args = reads.size() == 1 ? "-i ${reads[0]}" : "-i ${reads[0]} -I ${reads[1]}"
        def out_args = reads.size() == 1 ? "-o ${meta.output_dir}/fastq/${meta.id}.trimmed_r1.fastq.gz" 
                : "-o ${meta.output_dir}/fastq/${meta.id}.trimmed_r1.fastq.gz -O ${meta.output_dir}/fastq/${meta.id}.trimmed_r2.fastq.gz"
        def failed_args = meta.save_failed_trim ? ( 
            reads.size() == 1 ?
                "--failed_out ${meta.output_dir}/fastq/${meta.id}.failed.fastq.gz" 
                : "--failed_out ${meta.output_dir}/fastq/${meta.id}.failed.fastq.gz \
                    --unpaired1 ${meta.output_dir}/fastq/${meta.id}.unpaired_r1.fastq.gz \
                    --unpaired2 ${meta.output_dir}/fastq/${meta.id}.unpaired_r2.fastq.gz"
        ) : ""
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        // def output_dir = new File("${meta.id}/fastp")
        // if (!output_dir.exists()) {
        //     output_dir.mkdirs()
        // }
        """
            mkdir -p ${meta.output_dir}/fastq
            
            fastp \\
                ${args} \\
                -w ${task.cpus} \\
                ${input_args} \\
                ${out_args} \\
                -j ${meta.output_dir}/fastq/${meta.id}.json \\
                -h ${meta.output_dir}/fastq/${meta.id}.html \\
                ${failed_args}
        """
}