# Ref gen
ref_dir=/storage/home/duaghk/data/rna/reference/star
singularity exec \
-B /storage,/data \
/storage/images/star-2.7.11.sif \
STAR \
--runMode genomeGenerate \
--runThreadN 32 \
--genomeDir ${ref_dir} \
--genomeFastaFiles /storage/home/duaghk/data/rna/reference/Homo_sapiens_assembly38.basic.fasta \
--sjdbGTFfile /storage/home/duaghk/data/rna/reference/gencode.v47.basic.annotation.gtf \
--sjdbOverhang 100


# Onepass
singularity exec \
-B /storage,/data \
/storage/images/star-2.7.11.sif \
STAR \
--runMode alignReads \
--runThreadN 16 \
--genomeDir ${ref_dir} \
--outFileNamePrefix /storage/home/duaghk/data/rna/onepass/DU145. \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes All \
--outSAMattrRGline "ID:DU145 RG:DU145 SM:DU145_C0RTF PL:Illumina LB:SureSelect CN:CCLE" \
--readFilesIn \
/storage/home/duaghk/data/rna/fastq/DU_145_C0RTF_R1.fastq.gz \
/storage/home/duaghk/data/rna/fastq/DU_145_C0RTF_R2.fastq.gz \
--readFilesCommand zcat


# Genome generate for twopass
ref_dir=/storage/home/duaghk/data/rna/genome
singularity exec \
-B /storage,/data \
/storage/images/star-2.7.11.sif \
STAR \
--runMode genomeGenerate \
--runThreadN 32 \
--genomeDir ${ref_dir} \
--genomeFastaFiles /storage/home/duaghk/data/rna/reference/Homo_sapiens_assembly38.basic.fasta \
--sjdbGTFfile /storage/home/duaghk/data/rna/reference/gencode.v47.basic.annotation.gtf \
--sjdbOverhang 100 \
--sjdbFileChrStartEnd /storage/home/duaghk/data/rna/onepass/DU145.SJ.out.tab

# Twopass withinbam
ref_dir=/storage/home/duaghk/data/rna/genome
singularity exec \
-B /storage,/data \
/storage/images/star-2.7.11.sif \
STAR \
--runMode alignReads \
--runThreadN 16 \
--genomeDir ${ref_dir} \
--outFileNamePrefix /storage/home/duaghk/data/rna/twopass_withinbam/DU145. \
--outSAMattributes All \
--outSAMattrRGline "ID:DU145 RG:DU145 SM:DU145_C0RTF PL:Illumina LB:SureSelect CN:CCLE" \
--chimOutType WithinBAM \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn \
/storage/home/duaghk/data/rna/fastq/DU_145_C0RTF_R1.fastq.gz \
/storage/home/duaghk/data/rna/fastq/DU_145_C0RTF_R2.fastq.gz \
--readFilesCommand zcat

# RSEM prepare genome
ref_dir=/storage/home/duaghk/data/rna/reference/rsem/gencode_v47
singularity exec \
-B /storage,/data \
/storage/images/rsem-1.3.3.sif \
rsem-prepare-reference \
-p 16 \
--gtf /storage/home/duaghk/data/rna/reference/gencode.v47.basic.annotation.gtf \
/storage/home/duaghk/data/rna/reference/Homo_sapiens_assembly38.basic.fasta \
${ref_dir}

# RSEM quantification
singularity exec \
-B /storage,/data \
/storage/images/rsem-1.3.3.sif \
rsem-calculate-expression \
-p 16 \
--alignments \
--paired-end \
--strandedness reverse \
/storage/home/duaghk/data/rna/twopass_withinbam/DU145.Aligned.toTranscriptome.out.bam \
${ref_dir} \
/storage/home/duaghk/data/rna/rsem/DU145
