
nextflow.enable.dsl=2

include { ExtractInfoFromFastq } from '../../modules/utils'
include { Fastp } from '../../modules/fastp'
include { MEM } from '../../modules/bwa'
include { 
    samtools_Flagstat
    samtools_Sort
    samtools_Indexing
} from '../../modules/samtools'
include {
    gatk4_SortSam_aligned
    gatk4_SortSam_unaligned
    gatk4_FastqToSam
    gatk4_SamToFastq
    gatk4_MergeBamAlignment
    gatk4_CollectWgsMetrics
    gatk4_AddOrReplaceReadGroups
    gatk4_Mutect2
    gatk4_LearnReadOrientationModel
    gatk4_GetPileupSummaries
    gatk4_CalculateContamination
    gatk4_FilterMutectCalls
} from '../../modules/gatk4'
include {
    fgbio_ExtractUmisFromBam
    fgbio_GroupReadsByUmi
    fgbio_CallDuplexConsensusReads
    fgbio_FilterConsensusReads
    fgbio_CollectDuplexSeqMetrics
} from '../../modules/fgbio'

// For twice use
include { MEM as MEM2 } from '../../modules/bwa'
include {
    gatk4_SamToFastq as gatk4_SamToFastq2
    gatk4_MergeBamAlignment as gatk4_MergeBamAlignment2
} from '../../modules/gatk4'




workflow {
    // Define reference informations
    ch_ref = Channel.value(
            [
                fasta: "${params.ref_dir}/${params.ref_ver}/fasta/Homo_sapiens_assembly38.fasta",
                gnomad: "${params.ref_dir}/${params.ref_ver}/vcf/af-only-gnomad.hg38.vcf.gz",
                common_snp_intervals: "${params.ref_dir}/${params.ref_ver}/vcf/common_snp.bed",
                bait_bed: "/data/cfDNA/bed/Celemix/bait.bed",
                target_bed: "/data/cfDNA/bed/Celemix/target.bed",
            ]
    )
    // ch_ref.view()

    // Define default fastq informations
    ch_fastq = Channel
        .fromFilePairs(params.input_dir + '/*R{1,2}*.fastq.gz')
        .map {
            meta, fastq ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            // Set meta.single_end
            if (fastq.size() == 1) {
                fmeta.single_end = true
            } else {
                fmeta.single_end = false
            }
            [ fmeta, fastq ]
        }

    // Update Fastq channel's meta informations
    if (params.flowcell_id) {
        ch_fastq = ch_fastq | map { meta, reads -> 
                // Directly assign params.flowcell_id when it's provided
                meta.flowcell_id = params.flowcell_id
                meta.read_len = params.flowcell_id // Default value when ExtractInfoFromFastq is skipped
                [ meta, reads ]
            }
    } else {
        ch_fastq = ch_fastq | ExtractInfoFromFastq | map { meta, reads, v -> 
                // Extract flowcell_id and read_len from the ExtractInfoFromFastq output
                def (flowcell_id, read_len) = v.trim().split("\t")
                meta.flowcell_id = flowcell_id
                meta.read_len = read_len
                [ meta, reads ]
            }
    }

    ch_fastq = ch_fastq
        | map { meta, reads ->
            // Add other meta information
            meta.save_failed_trim = params.save_failed_trim
            meta.ref_ver = params.ref_ver
            meta.platform = params.platform
            meta.library = params.library
            meta.center = params.center
            // add UMI informations
            meta.umi_str = params.umi_str
            [ meta, reads ]
        }
    
    // Set module parameters.
    ch_fastq = ch_fastq
        | map { meta, reads ->
            // Add process arguments
            meta.args = [
                Fastp: """
                    --trim_poly_g \\
                    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \\
                    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
                """,
                MEM: "",
                // Samtools
                samtools_Flagstat: "-O json",
                // GATK4
                gatk4_FastqToSam: "",
                gatk4_MergeBamAlignment: """
                    -MAX_GAPS -1 \\
                    --CREATE_INDEX true \\
                    --CLIP_ADAPTERS false \\
                    -SO coordinate \\
                    --VALIDATION_STRINGENCY SILENT \\
                    --EXPECTED_ORIENTATIONS FR \\
                    --ALIGNER_PROPER_PAIR_FLAGS false
                """,
                gatk4_SamToFastq: "",
                gatk4_Mutect2: "-pairHMM AVX_LOGLESS_CACHING_OMP",
                gatk4_CollectWgsMetrics: "",
                // fgbio
                fgbio_ExtractUmisFromBam: """
                    --molecular-index-tags=ZA ZB \\
                    --single-tag=RX
                """,
                fgbio_GroupReadsByUmi: """
                    --raw-tag=RX \\
                    --min-map-q=10 \\
                    --edits=1
                """,
                fgbio_CallDuplexConsensusReads: """
                    --error-rate-pre-umi=45 \\
                    --error-rate-post-umi=30 \\
                    --min-input-base-quality=30 \\
                    --min-reads 2 1 1
                """,
                fgbio_FilterConsensusReads: """
                    --min-reads 2 1 1 \\
                    --max-base-error-rate 0.1 \\
                    --min-base-quality 20
                """,
                fgbio_CollectDuplexSeqMetrics: """
                    --duplex-umi-counts
                """,
            ]
            [ meta, reads ]
        }
    

    // Initialize output value
    // ch_fastq.view()

    // // workflow start
    // Primary: Make consensus reads
    // Set output dir
    ch_fastq = ch_fastq
        | map { meta, reads ->
            // Add output dir
            meta.output_dir = "${meta.id}/primary"
            meta.tag = "Primary"
            [ meta, reads ]
        }
    // println ch_fastq.view()

    // // trimming, Unaligned BAM
    ch_trimmed = ch_fastq | Fastp
    ch_unaligned = ch_trimmed.trimmed_reads | gatk4_FastqToSam | fgbio_ExtractUmisFromBam
    ch_umi_fastqs = ch_unaligned | gatk4_SamToFastq
    ch_aligned = MEM(ch_umi_fastqs, ch_ref)
    // ch_combined = ch_unaligned.join(ch_aligned, by: [0, 0])
    ch_combined = ch_unaligned.combine(ch_aligned)
        .filter { meta1, path1, meta2, path2 -> 
            meta1.id == meta2.id
        }
        .map { meta1, path1, meta2, path2 ->
            [meta1, path1, path2]
        }

    ch_merge = gatk4_MergeBamAlignment(ch_combined, ch_ref)
    samtools_Flagstat(ch_merge)
    ch_umi = ch_merge | fgbio_GroupReadsByUmi
    fgbio_CollectDuplexSeqMetrics(ch_umi)
    ch_consensus = ch_umi | fgbio_CallDuplexConsensusReads
    ch_filter = fgbio_FilterConsensusReads(ch_consensus, ch_ref)

    // // Second round
    // Set secondary output dir
    ch_filter = ch_filter
        | map { meta, reads ->
            // Add output dir
            meta.output_dir = "${meta.id}/secondary"
            meta.tag = "Secondary"
            [ meta, reads ]
        }
    // Update parameters.
    ch_filter = ch_filter
        | map { meta, reads ->
            meta.args = [
                MEM: "",
                // Samtools
                samtools_Flagstat: "-O json",
                // GATK4
                gatk4_SortSam_aligned: "-SO queryname",
                gatk4_SortSam_unaligned: "-SO queryname",
                gatk4_MergeBamAlignment: """
                    --VALIDATION_STRINGENCY=SILENT \\
                    --CLIP_ADAPTERS=false \\
                    --CREATE_INDEX=true \\
                    --EXPECTED_ORIENTATIONS=FR \\
                    --MAX_GAPS=-1 \\
                    --SORT_ORDER=coordinate \\
                    --ALIGNER_PROPER_PAIR_FLAGS=false \\
                    --ATTRIBUTES_TO_RETAIN=X0 \\
                    --ATTRIBUTES_TO_RETAIN=ZS \\
                    --ATTRIBUTES_TO_RETAIN=ZI \\
                    --ATTRIBUTES_TO_RETAIN=ZM \\
                    --ATTRIBUTES_TO_RETAIN=ZC \\
                    --ATTRIBUTES_TO_RETAIN=ZN \\
                    --ATTRIBUTES_TO_RETAIN=ad \\
                    --ATTRIBUTES_TO_RETAIN=bd \\
                    --ATTRIBUTES_TO_RETAIN=cd \\
                    --ATTRIBUTES_TO_RETAIN=ae \\
                    --ATTRIBUTES_TO_RETAIN=be \\
                    --ATTRIBUTES_TO_RETAIN=ce
                """,
                // SamToFastq2: """
                //     --INTERLEAVE true
                // """,
                gatk4_Mutect2: "-pairHMM AVX_LOGLESS_CACHING_OMP",
                gatk4_CollectWgsMetrics: "",
                gatk4_AddOrReplaceReadGroups: "--CREATE_INDEX true",
            ]
            [ meta, reads ]
        }
    ch_consensus_fastq = ch_filter | gatk4_SamToFastq2
    ch_consensus_unalign = ch_filter | gatk4_SortSam_unaligned
    ch_consensus_align = MEM2(ch_consensus_fastq, ch_ref)
    ch_consensus_align = ch_consensus_align | gatk4_SortSam_aligned
    ch_consensus = ch_consensus_unalign.combine(ch_consensus_align)
        .filter { meta1, path1, meta2, path2 -> 
            meta1.id == meta2.id
        }
        .map { meta1, path1, meta2, path2 ->
            [meta1, path1, path2]
        }

    ch_arbam = gatk4_MergeBamAlignment2(ch_consensus, ch_ref)
    ch_arbam = gatk4_AddOrReplaceReadGroups(ch_arbam)
    // Need to indexing using samtools ?
    gatk4_CollectWgsMetrics(ch_arbam, ch_ref)
    ch_pileup = gatk4_GetPileupSummaries(ch_arbam, ch_ref)
    (ch_cc, ch_segment) = gatk4_CalculateContamination(ch_pileup)
    (ch_mutect2, ch_f1r2) = gatk4_Mutect2(ch_arbam, ch_ref)
    ch_ob_prior = gatk4_LearnReadOrientationModel(ch_f1r2)
    // Combine all channelsabEF
    ch_mutect2_filter = ch_mutect2.combine(ch_ob_prior)
        .filter { meta1, mutect2_vcf, mutect2_vcf_idx, mutect2_vcf_stat, meta2, ob_prior -> 
            meta1.id == meta2.id
        }
        .map { meta1, mutect2_vcf, mutect2_vcf_idx, mutect2_vcf_stat, meta2, ob_prior -> 
            [meta1, mutect2_vcf, mutect2_vcf_idx, mutect2_vcf_stat, ob_prior]
        }
    ch_mutect2_filter = ch_mutect2_filter.combine(ch_segment)
        .filter { meta1, mutect2_vcf, mutect2_vcf_idx, mutect2_vcf_stat, ob_prior, meta2, segment_table -> 
            meta1.id == meta2.id            
        }
        .map { meta1, mutect2_vcf, mutect2_vcf_idx, mutect2_vcf_stat, ob_prior, meta2, segment_table -> 
            [meta1, mutect2_vcf, mutect2_vcf_idx, mutect2_vcf_stat, ob_prior, segment_table]
        }
    ch_mutect2_filter = ch_mutect2_filter.combine(ch_cc)
        .filter { meta1, mutect2_vcf, mutect2_vcf_idx, mutect2_vcf_stat, ob_prior, segment_table, meta2, cc_table -> 
            meta1.id == meta2.id            
        }
        .map { meta1, mutect2_vcf, mutect2_vcf_idx, mutect2_vcf_stat, ob_prior, segment_table, meta2, cc_table  -> 
            [meta1, mutect2_vcf, mutect2_vcf_idx, mutect2_vcf_stat, ob_prior, segment_table, cc_table]
        }
    ch_mutect2_filter = gatk4_FilterMutectCalls(ch_mutect2_filter, ch_ref)
    println "Mutect2 Filter"
    ch_mutect2_filter.view()
}


