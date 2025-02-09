nextflow.enable.dsl=2

process gatk4_FastqToSam {

    label 'middle_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(trimmed_reads)
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.unaligned.bam"), emit: unaligned_bam
    
    script:
        def rg_string = [
            meta.flowcell_id ? "-RG ${meta.flowcell_id}" : null,
            meta.id ? "-SM ${meta.id}" : null,
            meta.platform ? "-PL ${meta.platform}" : null,
            meta.library ? "-LB ${meta.library}" : null,
            meta.center ? "-CN ${meta.center}" : null,
        ].findAll { it }.join(" ")
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            FastqToSam \\
            ${args} \\
            ${rg_string} \\
            -F1 ${trimmed_reads[0]} \\
            -F2 ${trimmed_reads[1]} \\
            -O ${meta.output_dir}/bam/${meta.id}.unaligned.bam
        """
}

process gatk4_MergeBamAlignment {
    
    label 'middle_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(unaligned_bam), path(aligned_bam)
        val ref

    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.merged.bam"), emit: merged_bam
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            MergeBamAlignment \\
            ${args} \\
            -R ${ref.fasta} \\
            -UNMAPPED ${unaligned_bam} \\
            -ALIGNED ${aligned_bam} \\
            -O ${meta.output_dir}/bam/${meta.id}.merged.bam
        """
}

process gatk4_AddOrReplaceReadGroups {
    
    label 'standard'
    // tag "${meta.tag}"
   tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.RG-added.bam"), path("${meta.output_dir}/bam/${meta.id}.RG-added.bai"), emit: bam
    
    script:
        def rg_string = [
            meta.flowcell_id ? "--RGID ${meta.flowcell_id} --RGPU ${meta.flowcell_id}" : null,
            meta.id ? "--RGSM ${meta.id}" : null,
            meta.platform ? "--RGPL ${meta.platform}" : null,
            meta.library ? "--RGLB ${meta.library}" : null,
            meta.center ? "--RGCN ${meta.center}" : null,
        ].findAll { it }.join(" ")
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            AddOrReplaceReadGroups \\
            ${args} \\
            ${rg_string} \\
            -I ${bam} \\
            -O ${meta.output_dir}/bam/${meta.id}.RG-added.bam
        """
}

process gatk4_SamToFastq {
    
    label 'middle_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(unaligned_umi_bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/fastq/${meta.id}.UMI-tagged.R{1,2}.fastq.gz"), emit: umi_fastqs
    
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/fastq

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            SamToFastq \\
            ${args} \\
            -I ${unaligned_umi_bam} \\
            -F ${meta.output_dir}/fastq/${meta.id}.UMI-tagged.R1.fastq.gz \\
            -F2 ${meta.output_dir}/fastq/${meta.id}.UMI-tagged.R2.fastq.gz
        """
}

process gatk4_SortSam {
    
    label 'middle_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.${meta.prefix}.sorted.bam")
 
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            SortSam \\
            ${args} \\
            -I ${bam} \\
            -O ${meta.output_dir}/bam/${meta.id}.${meta.prefix}.sorted.bam
        """
}


process gatk4_SortSam_aligned {
    
    label 'middle_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.aligned.sorted.bam")
 
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            SortSam \\
            ${args} \\
            -I ${bam} \\
            -O ${meta.output_dir}/bam/${meta.id}.aligned.sorted.bam
        """
}

process gatk4_SortSam_unaligned {
    
    label 'middle_mem'
    // tag "${meta.tag}"
    tag 'Preprocessing'
    
    input: 
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${meta.output_dir}/bam/${meta.id}.unaligned.sorted.bam")
 
    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/bam

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            SortSam \\
            ${args} \\
            -I ${bam} \\
            -O ${meta.output_dir}/bam/${meta.id}.unaligned.sorted.bam
        """
}

// Metrics?
process gatk4_CollectWgsMetrics {

    label 'standard'
    tag 'Metrics'
    // tag "${meta.tag}"

    input: 
        tuple val(meta), path(bam), path(bai)
        val ref
    
    output:
        tuple val(meta), path("${meta.output_dir}/metric/${meta.id}.CollectWgsMetrics.txt"), emit: wgsmetric

    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/metric

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            CollectWgsMetrics \\
            ${args} \\
            -R ${ref.fasta} \\
            -I ${bam} \\
            -O ${meta.output_dir}/metric/${meta.id}.CollectWgsMetrics.txt
        """
}

// Mutect 2
process gatk4_Mutect2 {

    label 'middle_mem'
    // tag "${meta.tag}"
    tag 'VariantCall'

    input: 
        tuple val(meta), path(consensus_bam), path(consensus_bai)
        val ref
    
    output:
        tuple val(meta), path("${meta.output_dir}/variant/${meta.id}.mutect2.vcf.gz"), path("${meta.output_dir}/variant/${meta.id}.mutect2.vcf.gz.tbi"), path("${meta.output_dir}/variant/${meta.id}.mutect2.vcf.gz.stats"), emit: mutect2_vcf
        tuple val(meta), path("${meta.output_dir}/variant/${meta.id}.f1r2.tar.gz"), emit: f1r2

    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/variant

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            Mutect2 \\
            ${args} \\
            --native-pair-hmm-threads ${task.cpus} \\
            -R ${ref.fasta} \\
            --germline-resource ${ref.gnomad} \\
            -I ${consensus_bam} \\
            -O "${meta.output_dir}/variant/${meta.id}.mutect2.vcf.gz" \\
            --f1r2-tar-gz ${meta.output_dir}/variant/${meta.id}.f1r2.tar.gz
        """
}

process gatk4_LearnReadOrientationModel {

    label 'standard'
    tag 'VariantCall'
    // tag "${meta.tag}"

    input: 
        tuple val(meta), path(f1r2)
    
    output:
        tuple val(meta), path("${meta.output_dir}/variant/${meta.id}.artifact-prior.tar.gz"), emit: ob_prior

    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/variant

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            LearnReadOrientationModel \\
            ${args} \\
            -I ${f1r2} \\
            -O ${meta.output_dir}/variant/${meta.id}.artifact-prior.tar.gz \\
        """
}

process gatk4_GetPileupSummaries {

    label 'standard'
    tag 'VariantCall'
    // tag "${meta.tag}"

    input: 
        tuple val(meta), path(bam), path(bai)
        val ref
    
    output:
        tuple val(meta), path("${meta.output_dir}/variant/${meta.id}.pileups.table"), emit: pileup_summary

    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/variant

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            GetPileupSummaries \\
            ${args} \\
            -R ${ref.fasta} \\
            -I ${bam} \\
            -V ${ref.gnomad} \\
            -L ${ref.common_snp_intervals} \\
            -O ${meta.output_dir}/variant/${meta.id}.pileups.table
        """
}

process gatk4_CalculateContamination {

    label 'standard'
    tag 'VariantCall'
    // tag "${meta.tag}"

    input: 
        tuple val(meta), path(pileup_summary)
    
    output:
        tuple val(meta), path("${meta.output_dir}/variant/${meta.id}.cc.table"), emit: cc_table
        tuple val(meta), path("${meta.output_dir}/variant/${meta.id}.segments.table"), emit: segment_table

    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/variant

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            CalculateContamination \\
            ${args} \\
            -I ${pileup_summary} \\
            -O ${meta.output_dir}/variant/${meta.id}.cc.table \\
            -segments ${meta.output_dir}/variant/${meta.id}.segments.table
        """
}

process gatk4_FilterMutectCalls {

    label 'standard'
    tag 'VariantCall'
    // tag "${meta.tag}"

    input: 
        tuple val(meta), path(mutect2_vcf), path(mutect2_vcf_idx), path(mutect2_vcf_stat), path(ob_prior), path(segment_table), path(cc_table)
        val ref
    
    output:
        tuple val(meta), path("${meta.output_dir}/variant/${meta.id}.filtered.mutect2.vcf.gz"), emit: filtered_mutect2_vcf

    script:
        def args = meta.args.containsKey(task.process) ? meta.args[task.process].trim() : ""
        """
            mkdir -p ${meta.output_dir}/variant

            gatk \\
            --java-options \\
            \"-Xmx${task.memory.toGiga().toInteger()}g \\
            -XX:ConcGCThreads=${task.cpus} \\
            -Djava.io.tmpdir=${params.tmp_dir}\" \\
            FilterMutectCalls \\
            ${args} \\
            -R ${ref.fasta} \\
            -V ${mutect2_vcf} \\
            --stats ${mutect2_vcf_stat} \\
            -ob-priors ${ob_prior} \\
            --tumor-segmentation ${segment_table} \\
            --contamination-table ${cc_table} \\
            -O ${meta.output_dir}/variant/${meta.id}.filtered.mutect2.vcf.gz
        """
}


// process gatk4_MarkDuplicatesSpark {

//     label 'large_mem'
//     tag 'Preprocessing'
    
//     input: 
//         tuple val(meta), path(aligned_sam)
// }
