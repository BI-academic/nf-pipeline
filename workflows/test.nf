
process TestBed {
    output:
        stdout

    script:
        def base_args = """
                            --quantMode TranscriptomeSAM
                            --twopassMode Basic
                            --outSAMtype BAM Unsorted
                            --readFilesCommand zcat
                            --runRNGseed 0
                            --outFilterMultimapNmax 20
                            --alignSJDBoverhangMin 1
                            --outSAMattributes NH HI AS NM MD
                            --outSAMstrandField intronMotif
                        """.trim()
        """echo \"${base_args}\""""
}

process CreateDir {

    input:
        val meta

    output:
        path "${meta.id}/${meta.prefix}/${meta.id}"
    
    script:
        def args = meta.args[task.process]
        """            
            echo ${args} && mkdir -p ${meta.id}/${meta.prefix}/${meta.id}
        """        
}

workflow {
    // ch_input = Channel.value(
    //     [
    //         id: "UMI-Test", 
    //         prefix: "",
    //         args: [
    //             CreateDir: "ProcessName is CreateDir"
    //         ]
    //     ]
    // )
    // CreateDir(ch_input)


    channel_one = Channel.of(
        tuple(
            [
                id: "UMI_test1",
                unmet: true,
                prefix: "unaligned"
            ], 
            "unaligned.bam"
        )
    )
    // channel_one.view()

    channel_two = Channel.of(
        tuple(
            [
                id: "UMI_test1",
                unmet: true,
                prefix: "aligned"
            ], 
            "aligned.bam"
        )
    )
    // channel_two.view()
    // println channel_two

    // joined_ch = channel_one.join(channel_two, by: [0,0])
    joined_ch = channel_one.combine(channel_two)
        .filter { meta1, path1, meta2, path2 -> 
            meta1.id == meta2.id
        }
        .map { meta1, path1, meta2, path2 ->
            [meta1, path1, path2]
        }

    println joined_ch.view()

    // TestBed()
}