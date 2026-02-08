conda create -n bioinformatics python=3.9 -y
conda activate bioinformatics
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda hisat2 samtools sra-tools -y

mkdir exon_junction && cd exon_junction

prefetch SRR35927838
fastq-dump --split-files SRR35927838

mkdir -p hg38_index

cd hg38_index


curl -o hg38.tar.gz https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzvf hg38.tar.gz
cd ..
hisat2 -p 4 -x hg38_index/grch38/genome -U SRR35927838_1.fastq -S aligned_reads.sam
samtools sort -@ 4 -o cancer_mcf7.bam aligned_reads.sam
rm aligned_reads.sam
samtools index cancer_mcf7.bam
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
gunzip gencode.v45.annotation.gtf.gz

prefetch SRR25208720
fasterq-dump --split-3 SRR25208720 -e 4 -p
hisat2 -p 4 -x  hg38_index/grch38/genome -1 SRR25208720_1.fastq -2 SRR25208720_2.fastq -S normal_mcf10a.sam
samtools view -@ 4 -bS normal_mcf10a.sam | samtools sort -@ 4 -o normal_mcf10a_sorted.bam
samtools index normal_mcf10a_sorted.bam
rm normal_mcf10a.sam
