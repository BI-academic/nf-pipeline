nextflow run \
workflows/WTS \
--input_dir /storage/home/duaghk/data/rna/fastq \
--output_dir /storage/home/duaghk/data/rna/output \
--flowcell_id Test01 \
--read_len 101 \
--ref_dir /storage/home/duaghk/data/rna/reference \
--ref_ver hg38