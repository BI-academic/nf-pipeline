nextflow.enable.dsl=2

process ExtractInfoFromFastq {

    label 'standard'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path(reads), stdout

    script:
    """
    # Read the first two lines from the GZIP file using zhead
    first_line=\$(zcat ${reads[0]} | head -n 1)
    second_line=\$(zcat ${reads[0]} | head -n 2 | tail -n 1)

    # Extract the Flowcell ID from the first line
    flowcell_id=\$(echo "\$first_line" | awk -F ':' '{print \$3}')

    # Calculate the read length from the second line
    read_length=\$(echo -n "\$second_line" | wc -c)

    # Output the values
    echo "\${flowcell_id}\t\${read_length}"
    """
}


// workflow {
//     // Define test input
//     input_files = Channel
//         .fromPath('/data/fastq/TBD240401_20318_20240402/*.R1.fastq.gz')
//     input_files.view()
//     // Run the ExtractInfo process
//     input_files 
//         | ExtractInfoFromFastq 
//         | map { v -> 
//             def (flowcell_id, read_len) = v.trim().split("\t")
//             println "Flowcell ID: ${flowcell_id}, Read length: ${read_len}"
//             tuple(flowcell_id, read_len)
//         }
//         // | view()
// }