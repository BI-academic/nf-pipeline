
nextflow.enable.dsl=2

include { ExtractInfoFromFastq } from '../../modules/utils'

include { Fastp } from '../../modules/fastp'
include {
    STAR_OnePass
    STAR_TwoPass_withinbam
    STAR_TwoPass_chimeric
    STAR_GenomeGenerate
} from '../../modules/star'
include { RSEM_CalculateExpression } from '../../modules/rsem'


// Set Workflow
workflow {
    // Define reference channels
    ch_ref = Channel.value(
        [
            fasta: "${params.ref_dir}/${params.ref_ver}/Homo_sapiens_assembly38.basic.fasta",
            gtf: "${params.ref_dir}/${params.ref_ver}/gencode.v46.basic.annotation.gtf",
        ]
    )

    // Define Sample channels
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
            meta.ref_dir = params.ref_dir
            meta.platform = params.platform
            meta.library = params.library
            meta.center = params.center
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
                // STAR
                STAR_OnePass: """
                    --outSAMtype BAM SortedByCoordinate
                """,
                STAR_GenomeGenerate: "",
                STAR_TwoPass_withinbam: """
                    --chimOutType WithinBAM \\
                    --outSAMtype BAM SortedByCoordinate \\
                    --quantMode TranscriptomeSAM
                """,
                STAR_TwoPass_chimeric: """
                    --chimOutType Junctions \\
                    --outSAMtype BAM SortedByCoordinate
                """,
                // RSEM
                RSEM_CalculateExpression: """
                    --alignments \\
                    --paired-end \\
                    --strandedness reverse
                """
            ]
            [ meta, reads ]
        }

    // Set output directory
    ch_fastq = ch_fastq
        | map { meta, reads ->
            // Add output dir
            meta.output_dir = "${meta.id}"
            [ meta, reads ]
        }
    
    // Trimming first.
    ch_trimmed = ch_fastq | Fastp

    // Set STAR onepass output path
    ch_trimmed = ch_trimmed.trimmed_reads
        | map { meta, reads ->
            // Add output dir
            def genome_dir = "${meta.ref_dir}/${meta.ref_ver}/star"
            [ meta, reads, genome_dir ]
        }
    // Run Onepass
    ch_onepass = ch_trimmed | STAR_OnePass
    
    // Genome Generate
    ch_onepass = ch_onepass
        | map { meta, sj_out, bam, log_final, log, log_progress ->
            // replace genome dir to sample_specific
            def genome_dir = "${meta.id}/genome"
            [ meta, sj_out, genome_dir ]
        }

    ch_genome = ch_ref.combine(ch_onepass) | STAR_GenomeGenerate

    // Combine output
    ch_twopass_input = ch_trimmed.combine(ch_genome)
        .filter { meta1, reads, meta2, genome_dir -> 
            meta1.id == meta2.id
        }
        .map { meta1, reads, meta2, genome_dir -> 
            [meta1, reads, genome_dir]
        }

    ch_twopass_withinbam = ch_twopass_input | STAR_TwoPass_withinbam
    ch_twopass_chimeric = ch_twopass_input | STAR_TwoPass_chimeric
    ch_twopass_chimeric.view()

    ch_rsem = ch_twopass_withinbam
        | map { meta, transcript_bam, align_bam, sj_out, log_final, log, log_progress ->
            def rsem_ref = "${params.ref_dir}/${params.ref_ver}/rsem/gencode_v47"
            [meta, transcript_bam, rsem_ref]
        }
    ch_rsem = ch_rsem | RSEM_CalculateExpression
    ch_rsem.view()
}




