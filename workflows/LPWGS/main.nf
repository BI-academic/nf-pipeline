
nextflow.enable.dsl=2

include { ExtractInfoFromFastq } from '../../modules/utils'
include { Fastp } from '../../modules/fastp'
include { Dragmap } from '../../modules/dragmap'

workflow LPWGS {
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
    // ch_fastq.view()
    
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
            [ meta, reads ]
        }
        
    ch_fastq.view()

    // fastp_channel = ch_fastq | Fastp

    // // View the output for testing
    // dragmap_channel = fastp_channel.trimmed_reads | Dragmap

    // dragmap_channel.aligned_sam.view()
}

workflow {
    LPWGS()
}