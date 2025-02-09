
nextflow.enable.dsl=2

include { ExtractInfoFromFastq } from '../../modules/utils'

include {
    Align as OnePass
    Align as TwoPass
} from '../../modules/star'


// Set Workflow
workflow {
    // Define reference channels
    ch_ref = Channel.value(
        [
            fasta: "${params.ref_dir}/${params.ref_ver}/fasta/Homo_sapiens_assembly38.fasta",
            gtf: "${params.ref_dir}/${params.ref_ver}/gtf/gencode.v46.basic.annotation.gtf"
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
            meta.platform = params.platform
            meta.library = params.library
            meta.center = params.center
            // add UMI informations
            meta.umi_str = params.umi_str
            [ meta, reads ]
        }

}




