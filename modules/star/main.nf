nextflow.enable.dsl=2

process STAR_OnePass {

    label 'high'
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(trimmed_reads), path(genome_dir)
    
    output:
        tuple (
            val(meta), 
            path("${meta.id}/onepass/${meta.id}.SJ.out.tab"),
            path("${meta.id}/onepass/${meta.id}.Aligned.sortedByCoord.out.bam"), 
            path("${meta.id}/onepass/${meta.id}.Log.final.out"), 
            path("${meta.id}/onepass/${meta.id}.Log.out"), 
            path("${meta.id}/onepass/${meta.id}.Log.progress.out")
        )
        
    script:
        def rg_string = [
            meta.flowcell_id ? "ID:${meta.flowcell_id}" : null,
            meta.id ? "SM:${meta.id}" : null,
            meta.platform ? "PL:${meta.platform}" : null,
            meta.library ? "LB:${meta.library}" : null,
            meta.center ? "CN:${meta.center}" : null,
        ].findAll { it }.join(" ")
        
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        def read_command = trimmed_reads[0].toString().endsWith('.gz') ? '--readFilesCommand zcat' : ""
        """
            STAR \\
            --runMode alignReads \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${genome_dir} \\
            --outFileNamePrefix ${meta.id}/onepass/${meta.id}. \\
            --outSAMattrRGline ${rg_string} \\
            ${args} \\
            --readFilesIn ${trimmed_reads[0]} ${trimmed_reads[1]} \\
            ${read_command}
        """
}

process STAR_GenomeGenerate {

    label 'high'
    tag 'Preprocessing'
    
    input: 
        tuple path(fasta), path(gtf), val(meta), path(sjdb_path), path(genome_dir)

    output:
        tuple (
            val(meta), 
            path(genome_dir),
            path("${genome_dir}/Genome"),
            path("${genome_dir}/Log.out"),
            path("${genome_dir}/SA"),
            path("${genome_dir}/SAindex"),
            path("${genome_dir}/chrLength.txt"),
            path("${genome_dir}/chrName.txt"),
            path("${genome_dir}/chrNameLength.txt"),
            path("${genome_dir}/chrStart.txt"),
            path("${genome_dir}/exonGeTrInfo.tab"),
            path("${genome_dir}/exonInfo.tab"),
            path("${genome_dir}/geneInfo.tab"),
            path("${genome_dir}/genomeParameters.txt"),
            path("${genome_dir}/sjdbInfo.txt"),
            path("${genome_dir}/sjdbList.fromGTF.out.tab"),
            path("${genome_dir}/sjdbList.out.tab"),
            path("${genome_dir}/transcriptInfo.tab"),
        )
    
    script:
        def overhang = meta.read_len - 1
        def sjdb_args = sjdb_path ? "--sjdbFileChrStartEnd $sjdb_path": ''
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${genome_dir} \\
            --genomeFastaFiles ${fasta} \\
            --sjdbGTFfile ${gtf} \\
            --sjdbOverhang ${overhang} \\
            ${sjdb_args} \\
            ${args}
        """
}


process STAR_TwoPass_withinbam {

    label 'high'
    tag 'Processing'
    
    input: 
        tuple val(meta), path(trimmed_reads), path(genome_dir)
    
    output:
        tuple (
            val(meta), 
            path("${meta.id}/twopass_withinbam/${meta.id}.Aligned.toTranscriptome.out.bam"), 
            path("${meta.id}/twopass_withinbam/${meta.id}.Aligned.sortedByCoord.out.bam"), 
            path("${meta.id}/twopass_withinbam/${meta.id}.SJ.out.tab"), 
            path("${meta.id}/twopass_withinbam/${meta.id}.Log.final.out"), 
            path("${meta.id}/twopass_withinbam/${meta.id}.Log.out"), 
            path("${meta.id}/twopass_withinbam/${meta.id}.Log.progress.out")
        )
        
    script:
        def rg_string = [
            meta.flowcell_id ? "ID:${meta.flowcell_id}" : null,
            meta.id ? "SM:${meta.id}" : null,
            meta.platform ? "PL:${meta.platform}" : null,
            meta.library ? "LB:${meta.library}" : null,
            meta.center ? "CN:${meta.center}" : null,
        ].findAll { it }.join(" ")
        
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        def read_command = trimmed_reads[0].toString().endsWith('.gz') ? '--readFilesCommand zcat' : ""
        """
            STAR \\
            --runMode alignReads \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${genome_dir} \\
            --outFileNamePrefix ${meta.id}/twopass_withinbam/${meta.id}. \\
            --outSAMattrRGline ${rg_string} \\
            ${args} \\
            --readFilesIn ${trimmed_reads[0]} ${trimmed_reads[1]} \\
            ${read_command}
        """
}

process STAR_TwoPass_chimeric {

    label 'high'
    tag 'Processing'
    
    input: 
        tuple val(meta), path(trimmed_reads), path(genome_dir)
    
    output:
        tuple (
            val(meta), 
            path("${meta.id}/twopass_chimeric/${meta.id}.Aligned.sortedByCoord.out.bam"), 
            path("${meta.id}/twopass_chimeric/${meta.id}.SJ.out.tab"), 
            path("${meta.id}/twopass_chimeric/${meta.id}.Log.final.out"), 
            path("${meta.id}/twopass_chimeric/${meta.id}.Log.out"), 
            path("${meta.id}/twopass_chimeric/${meta.id}.Log.progress.out")
        )
        
    script:
        def rg_string = [
            meta.flowcell_id ? "ID:${meta.flowcell_id}" : null,
            meta.id ? "SM:${meta.id}" : null,
            meta.platform ? "PL:${meta.platform}" : null,
            meta.library ? "LB:${meta.library}" : null,
            meta.center ? "CN:${meta.center}" : null,
        ].findAll { it }.join(" ")
        
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        def read_command = trimmed_reads[0].toString().endsWith('.gz') ? '--readFilesCommand zcat' : ""
        """
            STAR \\
            --runMode alignReads \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${genome_dir} \\
            --outFileNamePrefix ${meta.id}/twopass_chimeric/${meta.id}. \\
            --outSAMattrRGline ${rg_string} \\
            ${args} \\
            --readFilesIn ${trimmed_reads[0]} ${trimmed_reads[1]} \\
            ${read_command}
        """
}


