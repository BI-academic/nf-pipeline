nextflow.enable.dsl=2

process Align {

    label 'high'
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(trimmed_reads)
        val ref
    
    output:
        tuple val(meta), path('*Log.final.out'), emit: log_final
        tuple val(meta), path('*Log.out'), emit: log_out
        tuple val(meta), path('*Log.progress.out'), emit: log_progress

        tuple val(meta), path('*d.out.bam'), optional:true, emit: bam
        tuple val(meta), path("${prefix}.sortedByCoord.out.bam"), optional:true, emit: bam_sorted
        tuple val(meta), path("${prefix}.Aligned.sortedByCoord.out.bam"), optional:true, emit: bam_sorted_aligned
        tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
        tuple val(meta), path('*Aligned.unsort.out.bam'), optional:true, emit: bam_unsorted
        tuple val(meta), path('*fastq.gz'), optional:true, emit: fastq
        tuple val(meta), path('*.tab'), optional:true, emit: tab
        tuple val(meta), path('*.SJ.out.tab'), optional:true, emit: spl_junc_tab
        tuple val(meta), path('*.ReadsPerGene.out.tab'), optional:true, emit: read_per_gene_tab
        tuple val(meta), path('*.out.junction'), optional:true, emit: junction
        tuple val(meta), path('*.out.sam'), optional:true, emit: sam
        tuple val(meta), path('*.wig'), optional:true, emit: wig
        tuple val(meta), path('*.bg'), optional:true, emit: bedgraph

    script:
        def rg_string = [
            meta.flowcell_id ? "RG:${meta.flowcell_id}" : null,
            meta.id ? "SM:${meta.id}" : null,
            meta.platform ? "PL:${meta.platform}" : null,
            meta.library ? "LB:${meta.library}" : null,
            meta.center ? "CN:${meta.center}" : null,
        ].findAll { it }.join(" ")
        
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        def read_command = trimmed_reads[0].toString().endsWith('.gz') ? '--readFilesCommand zcat -' : ""
        """
            STAR \\
            --runMode alignReads \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${ref.genome_dir} \\
            --outFileNamePrefix ${meta.id}/onepass/${meta.id}. \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMattrRGline ${rg_string} \\
            --readFilesIn ${trimmed_reads[0]} ${trimmed_reads[1]} \\
            ${read_command}
        """
}

process GenomeGenerate {

    label 'high'
    tag 'Reference'
    
    input: 
        tuple path(fasta), path(gtf), val(read_len), path(ref_dir)
        path sjdb_path
        path twopass_genome

    output:
        path genome_dir
    
    def filter = opt.name != 'NO_FILE' ? "--filter $opt" : ''
    script:
        def overhang = read_len - 1
        def genome_dir =  twopass_genome ? twopass_genome : ref_dir
        def sjdb_args = sjdb_path ? "--sjdbFileChrStartEnd $sjdb_path ": ''
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${genome_dir} \\
            --genomeFastaFiles ${fasta} \\
            --sjdbGTFfile ${gtf} \\
            --sjdbOverhang ${overhang} \\
            ${sdjb_args} \\
            ${args}
        """
}

