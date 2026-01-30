#!/bin/bash
###The complete analytical pipeline, including software tools and parameter configurations, is provided here.####

###Identify one-to-one orthologous genes. The longest protein sequence from each species was chosen for downstream analyses.
python orthofinder.py -f ./path -S diamond -t 5 -a 5

###Identify BUSCO genes.
busco -i pep.fa -l metazoa_odb10 -o ./path -m proteins

##ω.
###Identify one-to-one orthologous genes for group species.
python orthofinder.py -f path -S diamond -t 5 -a 5
/home/muscle -align -in pep.fa -output out.fa
pal2nal.pl out.fa cds.fa -output fasta -nogap -nomismatch > pal2nal_out.fa
###concatenated the alignments within each group.
codeml file.ctl ###For detailed parameters, see the Methods section.

###STAR index
STAR --runThreadN 60 --genomeSAindexNbases 10 --runMode genomeGenerate --genomeDir index1 --genomeFastaFiles genomic.fa --sjdbGTFfile genomic.gtf -sjdbOverhang 99

###bowtie2 index
bowtie2-build genomic.fa ./bowtie2_index/index


####RNA-seq data analysis pipeline
###For paired layout
fastp --thread 16 -i ${i}_1.fastq.gz -o ./path/${i}_1.fastq.gz -I ${i}_2.fastq.gz -O ./path/${i}_2.fastq.gz
###For single layout
fastp --thread 16 -i ${i}.fastq.gz -o ./path/${i}.fastq.gz
###STAR mapping
STAR --runThreadN 10 --genomeDir STAR_index --readFilesCommand zcat --outSAMunmapped None --outFilterMultimapNmax 1 --alignEndsType EndToEnd --outSAMstrandField intronMotif --outSJfilterReads Unique --readFilesIn ${i}.fastq.gz --outFileNamePrefix ./${i}
###sort sam file
samtools view -@ 10 -S ${i}Aligned.out.sam -b > ${i}.bam
samtools sort -@ 10 ${i}.bam -o ${i}_sorted.bam
samtools index -@ 10 ${i}_sorted.bam

###Assembly RNA-seq gtf file
stringtie -p 10 bam_file -G stringtie.gtf -o RNAseq.gtf
####The stringtie.gtf file provides the start and end coordinates of genes and transcripts based on the downloaded genome annotation file.
gtf_to_alignment_gff3.pl RNAseqRef.gtf > RNAseqRef.gff3
gtf_genome_to_cdna_fasta.pl RNAseqRef.gtf genomic.fa > RNAseqRef.fasta
TransDecoder.LongOrfs -t RNAseqRef.fasta
TransDecoder.Predict -t RNAseqRef.fasta
cdna_alignment_orf_to_genome_orf.pl RNAseqRef.fasta.transdecoder.gff3 RNAseqRef.gff3 RNAseqRef.fasta > RNAseqRefGenome.gff3

###gene counts
featureCounts -T 64 -p -t gene -g gene_id -a genomic.gtf -o read_counts.txt bam_file


###ATI pipeline
###SEASTAR (python 2.7): Systematic Evaluation of Alternative STArt site in RNA.
bash SEASTAR.sh -A bam_file -B bam_file -o path -g RNAseqRefGenome.gff3 -i sorted_reference.genome.size.txt -s bowtie2_index/index.fa -p 5 -d 30
###The parameter -g can be used to run SEASTAR in a De novo mode.
###Based on the FilteredNrtss.annotation file, we computed diversity using an R script; details of the calculation are provided in the Methods section.


###APA pipeline
###TAPAS: tool for alternative polyadenylation site analysis. Bioinformatics, 2018.
###Command to calculate coverage from bam file:
samtools depth bam_file > read_coverage.txt
gtfToGenePred -genePredExt RNAseqRefGenome.gff3 -ignoreGroupsWithoutExons refFlat_sf.txt
./TAPAS_master/Finding_APA_Sites/APA_sites_detection -ref refFlat_sf.txt -cov read_coverage.txt -l read_length -o expression_with_read.txt
###Based on the expression_with_read.txt file, we computed diversity using an R script; details of the calculation are provided in the Methods section.


###AS pipeline
###Based on the SJ.out.tab file, we computed diversity using an R script; details of the calculation are provided in the Methods section.


###CAGE-seq pipeline.
###Extract all reads starting with G.
seqkit grep -s -r -p "^G" fastq.gz | gzip > fastq.gz
###hisat2 index.
hisat2-build -p 16 genomic.fa hisat2Index
###Remove rRNA-related reads. rRNA.fa from ENSEMBL genome file.
rRNAdust -t 60 -e 3 rRNA.fa fastq > .fastq
###mapping.
hisat2 -x hisat2Index -p 60 --no-softclip -U fastq -S sam
###Unique mapping reads.
samtools view -@ 10 -bS -q 20 sam > unique.sam
samtools sort -@ 10 unique.bam -o sorted.bam


###3'-end-seq pipeline.
###bowtie2 mapping
bowtie2 -p 60 -x index -1 _1.fastq.gz -2 _2.fastq.gz -S sam -L 25 -N 0 -i S,1,1.15 --local -5 4
###APA site identification was subsequently carried out using an in-house analysis script.

###All other related analyses were conducted using custom R or Python scripts developed in-house.