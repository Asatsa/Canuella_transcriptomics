#!/bin/sh

## path to raw data
/home/anabwire/data/rawdata/A006850263_188710_S55_L000_R2_001.fastq.gz
/home/anabwire/data/rawdata/A006850263_188705_S53_L000_R1_001.fastq.gz
/home/anabwire/data/rawdata/A006850263_188722_S61_L000_R2_001.fastq.gz
/home/anabwire/data/rawdata/A006850263_188718_S59_L000_R2_001.fastq.gz
/home/anabwire/data/rawdata/A006850263_188716_S58_L000_R2_001.fastq.gz
/home/anabwire/data/rawdata/A006850263_188714_S57_L000_R2_001.fastq.gz
/home/anabwire/data/rawdata/A006850263_188724_S62_L000_R2_001.fastq.gz
/home/anabwire/data/rawdata/A006850263_188720_S60_L000_R1_001.fastq.gz
/home/anabwire/data/rawdata/A006850263_188708_S54_L000_R2_001.fastq.gz
/home/anabwire/data/rawdata/A006850263_188712_S56_L000_R2_001.fastq.gz

#Quality Control
##rcorrector on raw data
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188710_S55_R1.fastq.gz,/home/anabwire/data/raw_data/188712_S56_R1.fastq.gz,/home/anabwire/data/raw_data/188714_S57_R1.fastq.gz,/home/anabwire/data/raw_data/188716_S58_R1.fastq.gz,/home/anabwire/data/raw_data/188718_S59_R1.fastq.gz,/home/anabwire/data/raw_data/188720_S60_R1.fastq.gz,/home/anabwire/data/raw_data/188722_S61_R1.fastq.gz,/home/anabwire/data/raw_data/188724_S62_R1.fastq.gz, -2 /home/anabwire/data/raw_data/188710_S55_R2.fastq.gz,/home/anabwire/data/raw_data/188712_S56_R2.fastq.gz,/home/anabwire/data/raw_data/188714_S57_R2.fastq.gz,/home/anabwire/data/raw_data/188716_S58_R2.fastq.gz,/home/anabwire/data/raw_data/188718_S59_R2.fastq.gz,/home/anabwire/data/raw_data/188720_S60_R2.fastq.gz,/home/anabwire/data/raw_data/188724_S62_R2.fastq.gz,/home/anabwire/data/raw_data/188722_S61_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrectorall.out &

##rcorrector on trimmed data
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/trimmed_reads/A006850263_188705_S53_L000_R1_001_val_1.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188708_S54_L000_R1_001_val_1.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188710_S55_L000_R1_001_val_1.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188712_S56_L000_R1_001_val_1.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188714_S57_L000_R1_001_val_1.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188716_S58_L000_R1_001_val_1.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188718_S59_L000_R1_001_val_1.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188720_S60_L000_R1_001_val_1.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188722_S61_L000_R1_001_val_1.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188724_S62_L000_R1_001_val_1.fq.gz -2 /home/anabwire/data/trimmed_reads/A006850263_188705_S53_L000_R2_001_val_2.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188708_S54_L000_R2_001_val_2.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188710_S55_L000_R2_001_val_2.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188712_S56_L000_R2_001_val_2.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188714_S57_L000_R2_001_val_2.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188716_S58_L000_R2_001_val_2.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188718_S59_L000_R2_001_val_2.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188720_S60_L000_R1_001_val_2.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188722_S61_L000_R2_001_val_2.fq.gz,/home/anabwire/data/trimmed_reads/A006850263_188724_S62_L000_R2_001_val_2.fq.gz -od /home/anabwire/data/corrected_reads -t 50 >> rcorrectorall1.out &

##rcorrector on rawdata individual file pairs
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188705_S53_R1.fastq.gz -2 /home/anabwire/data/raw_data/188705_S53_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 16 >> rcorrector.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188708_S54_R1.fastq.gz -2 /home/anabwire/data/raw_data/188708_S54_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrector1.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188710_S55_R1.fastq.gz -2 /home/anabwire/data/raw_data/188710_S55_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrector55.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188712_S56_R1.fastq.gz -2 /home/anabwire/data/raw_data/188712_S56_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrector56.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188714_S57_R1.fastq.gz -2 /home/anabwire/data/raw_data/188714_S57_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrector57.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188716_S58_R1.fastq.gz -2 /home/anabwire/data/raw_data/188716_S58_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrector58.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188718_S59_R1.fastq.gz -2 /home/anabwire/data/raw_data/188718_S59_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrector59.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188720_S60_R1.fastq.gz -2 /home/anabwire/data/raw_data/188720_S60_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrector60.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188722_S61_R1.fastq.gz -2 /home/anabwire/data/raw_data/188722_S61_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrector61.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/raw_data/188724_S62_R1.fastq.gz -2 /home/anabwire/data/raw_data/188724_S62_R2.fastq.gz -od /home/anabwire/data/corrected_reads -t 20 >> rcorrector62.out &
nohup perl /home/anabwire/rcorrector/run_rcorrector.pl -1 /home/anabwire/data/trimmed_reads/A006850263_188720_S60_L000_R1_001_val_1.fq.gz -2 /home/anabwire/data/trimmed_reads/A006850263_188720_S60_L000_R2_001_val_2.fq.gz -od /home/anabwire/data/corrected_reads -t 50 >> rcorrector60.out &


##Quality control using trim galore a dependence of Cut adapt


## trim galore on corrected reads
nohup trim_galore -a "file:/home/anabwire/data/raw_data/adapters_read1.fasta" -a2 "file:/home/anabwire/data/raw_data/adapters_read2.fasta" --paired --retain_unpaired --phred33 --output_dir trimmed_reads -q 15 --stringency 5 -e 0.05 --length 36 --fastqc 188705_S53_R1.cor.fq.gz 188705_S53_R2.cor.fq.gz 188708_S54_R1.cor.fq.gz 188708_S54_R2.cor.fq.gz 188710_S55_R1.cor.fq.gz 188710_S55_R2.cor.fq.gz 188712_S56_R1.cor.fq.gz 188712_S56_R2.cor.fq.gz 188714_S57_R1.cor.fq.gz 188714_S57_R2.cor.fq.gz 188716_S58_R1.cor.fq.gz 188716_S58_R2.cor.fq.gz 188718_S59_R1.cor.fq.gz 188718_S59_R2.cor.fq.gz 188720_S60_R1.cor.fq.gz 188720_S60_R2.cor.fq.gz 188722_S61_R1.cor.fq.gz 188722_S61_R2.cor.fq.gz 188724_S62_R1.cor.fq.gz 188724_S62_R2.cor.fq.gz >> rcorrectortrim.out &

## trim galore on raw data
nohup trim_galore -a "file:/home/anabwire/data/raw_data/adapters_read1.fasta" -a2 "file:/home/anabwire/data/raw_data/adapters_read2.fasta" --paired --retain_unpaired --phred33 --output_dir trimmed_reads -q 15 --stringency 5 -e 0.05 --length 36 --fastqc /home/anabwire/data/rawdata/A006850263_188705_S53_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188705_S53_L000_R2_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188708_S54_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188708_S54_L000_R2_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188710_S55_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188710_S55_L000_R2_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188712_S56_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188712_S56_L000_R2_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188714_S57_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188714_S57_L000_R2_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188716_S58_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188716_S58_L000_R2_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188718_S59_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188718_S59_L000_R2_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188720_S60_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188720_S60_L000_R2_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188722_S61_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188722_S61_L000_R2_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188724_S62_L000_R1_001.fastq.gz /home/anabwire/data/rawdata/A006850263_188724_S62_L000_R2_001.fastq.gz >> fu1.out &

## To be clear the sequence of events for was rcorrector then trim galore then trinity assembly



### trinity
 Trinity --seqType fq --left reads_1.fq --right reads_2.fq --CPU 6 --max_memory 20G ## the basic command


## used absolute file paths
                nohup Trinity --seqType fq \
                        --left /home/anabwire/data/corrected_reads/A006850263_188705_S53_L000_R1_001_val_1.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188708_S54_L000_R1_001_val_1.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188710_S55_L000_R1_001_val_1.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188712_S56_L000_R1_001_val_1.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188714_S57_L000_R1_001_val_1.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188716_S58_L000_R1_001_val_1.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188718_S59_L000_R1_001_val_1.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188720_S60_L000_R1_001_val_1.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188722_S61_L000_R1_001_val_1.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188724_S62_L000_R1_001_val_1.cor.fq.gz \
                        --right /home/anabwire/data/corrected_reads/A006850263_188705_S53_L000_R2_001_val_2.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188708_S54_L000_R2_001_val_2.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188710_S55_L000_R2_001_val_2.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188712_S56_L000_R2_001_val_2.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188714_S57_L000_R2_001_val_2.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188716_S58_L000_R2_001_val_2.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188718_S59_L000_R2_001_val_2.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188720_S60_L000_R2_001_val_2.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188722_S61_L000_R2_001_val_2.cor.fq.gz,/home/anabwire/data/corrected_reads/A006850263_188724_S62_L000_R2_001_val_2.cor.fq.gz \
                        --CPU 15 --max_memory 85G >> trinity6.out &
##Installing tools
                JELLYFISH_VERSION 2.3.0
                RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz && \
                    tar xvf jellyfish-2.3.0.tar.gz && \
                    cd jellyfish-2.3.0/ && \
                    ./configure && make && make install

                    SALMON_VERSION=1.5.2
                    RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/Salmon-1.5.2_linux_x86_64.tar.gz && \
                        tar xvf Salmon-1.5.2_linux_x86_64.tar.gz && \
                        ln -s $SRC/salmon-1.5.2_linux_x86_64/bin/salmon $BIN/.

### Trinity assembly quality assessment via bowtie

## build a bowtie2 index for the Trinity assembly, required before running the alignment:
noshup bowtie2-build trinity_out_dir/Trinity.fasta trinity_out_dir/Trinity.fasta >> bowtie.out&


## align the reads to the assembly

### used individual file pairs
      bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 reads_1.fq -2 reads_2.fq  \
     2>align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam

    nohup bowtie2 -p 10 -q --no-unal -k 20 -x ../trinity_out_dir/Trinity.fasta -1 ../trimmed_reads/A006850263_188705_S53_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188705_S53_L000_R2_001_val_2.fq.gz  \
     2>align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam >> bowtie53.out &
     nohup bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 ../trimmed_reads/A006850263_188705_S53_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188705_S53_L000_R2_001_val_2.fq.gz  \
     2>align_stats1.txt| samtools view -@10 -Sb -o bowtie2.bam >> bowtie53.out &

     nohup bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 ../trimmed_reads/A006850263_188708_S54_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188708_S54_L000_R2_001_val_2.fq.gz  \
     2>align_stats2.txt| samtools view -@10 -Sb -o bowtie2.bam >> bowtie54.out &

     nohup bowtie2 -p 10 -q --no-unal -k 20 -x ../trinity_out_dir/Trinity.fasta -1 ../trimmed_reads/A006850263_188710_S55_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188710_S55_L000_R2_001_val_2.fq.gz  \
     2>align_stats3.txt| samtools view -@10 -Sb -o bowtie3.bam >> bowtie55.out &

     nohup bowtie2 -p 10 -q --no-unal -k 20 -x ../../trinity_out_dir/Trinity.fasta -1 ../../trimmed_reads/A006850263_188712_S56_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188712_S56_L000_R2_001_val_2.fq.gz  \
     2>align_stats4.txt| samtools view -@10 -Sb -o bowtie4.bam >> bowtie56.out &

     nohup bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 ../trimmed_reads/A006850263_188714_S57_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188714_S57_L000_R2_001_val_2.fq.gz  \
     2>align_stats5.txt| samtools view -@10 -Sb -o bowtie5.bam >> bowtie57.out &

     nohup bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 ../trimmed_reads/A006850263_188716_S58_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188716_S58_L000_R2_001_val_2.fq.gz  \
     2>align_stats6.txt| samtools view -@10 -Sb -o bowtie6.bam >> bowtie58.out &

     nohup bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 ../trimmed_reads/A006850263_188718_S59_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188718_S59_L000_R2_001_val_2.fq.gz  \
     2>align_stats7.txt| samtools view -@10 -Sb -o bowtie7.bam >> bowtie59.out &

     nohup bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 ../trimmed_reads/A006850263_188720_S60_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188720_S60_L000_R2_001_val_2.fq.gz  \
     2>align_stats8.txt| samtools view -@10 -Sb -o bowtie8.bam >> bowtie60.out &

     nohup bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 ../trimmed_reads/A006850263_188722_S61_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188722_S61_L000_R2_001_val_2.fq.gz  \
     2>align_stats9.txt| samtools view -@10 -Sb -o bowtie9.bam >> bowtie61.out &

     nohup bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 ../trimmed_reads/A006850263_188724_S62_L000_R1_001_val_1.fq.gz -2 ../trimmed_reads/A006850263_188724_S62_L000_R2_001_val_2.fq.gz  \
     2>align_stats10.txt| samtools view -@10 -Sb -o bowtie10.bam >> bowtie62.out &

     samtools sort bowtie2.bam -o bowtie2.coordSorted.bam
     samtools index bowtie2.coordSorted.bam
     samtools faidx Trinity.fasta
     igv.sh -g Trinity.fasta  bowtie2.coordSorted.bam



#     Busco analysis

$ conda create -n env2 busco
$ source activate env2
generate_busco_config.py > config.ini
BUSCO_CONFIG_FILE=config.ini run_busco



nohup busco -m transcriptome -i Trinity.fasta -o busco_arthropoda -l arthropoda_odb10 >> busc.out &


  #check
  cat run_nema_busco_metazoa/short_summary_nema_busco_metazoa.txt

## Transcript abundance estimation
## Just prepare the reference for alignment and abundance estimation
/home/anabwire/miniconda3/envs/trinity/bin/align_and_estimate_abundance.pl --transcripts Trinity.fasta --samples_file Samples.tsv --est_method salmon --trinity_mode --prep_reference

## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)

/home/anabwire/miniconda3/envs/trinity/bin/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --output_dir rsem_outdir
#
## prep the reference and run the alignment/estimation

#### i USED this code
/home/anabwire/miniconda3/envs/trinity/bin/align_and_estimate_abundance.pl --seqType fq  \
    --samples_file samples.txt  --transcripts Trinity.fasta \
    --est_method salmon  --trinity_mode   --prep_reference
###https://bioinformaticsdotca.github.io/rnaseq_2018_tutorial6

### it worked this is from the trinity wiki pAGE
/home/anabwire/miniconda3/envs/trinity/bin/abundance_estimates_to_matrix.pl --est_method salmon \
      --gene_trans_map Trinity.fasta.gene_trans_map \
      --quant_files quant_files.list \
      --name_sample_by_basedir

### did not work
/home/anabwire/miniconda3/envs/trinity/bin/count_matrix_features_given_MIN_TPM_threshold.pl \
          genes_matrix.TPM.not_cross_norm | tee genes_matrix.TPM.not_cross_norm.counts_by_min_TPM


#and

/home/anabwire/miniconda3/envs/trinity/bin/count_matrix_features_given_MIN_TPM_threshold.pl \
          trans_matrix.TPM.not_cross_norm | tee trans_matrix.TPM.not_cross_norm.counts_by_min_TPM
###Contig Ex90N50 Statistic and Ex90 Gene Count
/home/anabwire/miniconda3/envs/trinity/bin/contig_ExN50_statistic.pl \
     Trinity.isoform.TMM.EXPR.matrix Trinity.fasta transcript | tee ExN50.transcript.stats

/home/anabwire/miniconda3/envs/trinity/bin/plot_ExN50_statistic.Rscript  ExN50.transcript.stats

     xpdf ExN50.transcript.stats.plot.pdf ## a code for viewing the pdf plot from the commandline but it didnt waork for me

### how to tell how many corresponds to the E90 peak
cat Trinity.isoform.TMM.EXPR.matrix.E-inputs|  egrep -v ^\# | awk '$1 <= 90' | wc -l

###OTHER ExN50PLOTS
/home/anabwire/miniconda3/envs/trinity/bin/contig_ExN50_statistic.pl \
        Trinity.isoform.TMM.EXPR.matrix Trinity.fasta gene | tee ExN50.gene.stats

#followed by plotting:

/home/anabwire/miniconda3/envs/trinity/bin/plot_ExN50_statistic.Rscript  ExN50.gene.stats

   xpdf ExN50.gene.stats.plot.pdf

##Estimating TPM thresholds for transcript counting and filtering i.e Ex90

/home/anabwire/miniconda3/envs/trinity/bin/try_estimate_TPM_filtering_threshold.Rscript --E_inputs Trinity.isoform.TMM.EXPR.matrix.E-inputs

    xpdf estimate_TPM_threshold.pdf

# this step is not working for some reason, but its just for estimating count of the the minimum expression levels.

## this is a code for downloading the pdf plots I have generated
rsync -PauL anabwire@we11sv01:/home/anabwire/data/aquality/ExN50.transcript.stats.plot.pdf /Users/nabwireasatsa/Thesis/Canuella_transcriptomics/results

##Quality Check Your Samples and Biological Replicates
### Compare replicates for each of your samples
$TRINITY_HOME/Analysis/DifferentialExpression/PtR --matrix counts.matrix \
                --samples samples.txt --log2 --CPM \
                --min_rowSums 10 \
                --compare_replicates
###Compare Replicates Across Samples
                $TRINITY_HOME/Analysis/DifferentialExpression/PtR \
                          --matrix Trinity_trans.counts.matrix \
                          --min_rowSums 10 \
                          -s samples.txt --log2 --CPM --sample_cor_matrix
### PCA plot
$TRINITY_HOME/Analysis/DifferentialExpression/PtR \
    --matrix Trinity_trans.counts.matrix \
    -s samples.txt --min_rowSums 10 --log2 \
    --CPM --center_rows \
    --prin_comp 3

###DE ANALYSIS FROM HERE WE PURELY USE THE TURORIAL
nohup /home/jboyen/bin/trinityrnaseq-v2.15.1/Analysis/DifferentialExpression/run_DE_analysis.pl \
      --matrix ../Trinity.isoform.counts.matrix \
      --samples_file samples.txt \
      --method DESeq2 \
      --output DESeq2_trans >> DE.out &

##DE using DESeq2 did not work in trinity, it was done separately on local R


/usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl

      $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
            --matrix Trinity.isoform.counts.matrix \
            --samples_file samples.txt \
            --method DESeq2 \
            --output DESeq2_trans
#DE FOR TRANSCRIPTS
/home/anabwire/miniconda3/envs/trinity/bin/run_DE_analysis.pl \
                      --matrix ../Trinity.isoform.counts.matrix \
                      --method edgeR \
                      --samples_file ../samples.txt
/home/anabwire/miniconda3/envs/trinity/bin/analyze_diff_expr.pl --matrix ../../Trinity.isoform.TMM.EXPR.matrix --samples ../samples.txt -P 1e-3 -C 2
## DE for genes
/home/anabwire/miniconda3/envs/trinity/bin/run_DE_analysis.pl \
                        --matrix ../Trinity.gene.counts.matrix \
                        --method edgeR \
                        --samples_file ../samples.txt

/home/anabwire/miniconda3/envs/trinity/bin/analyze_diff_expr.pl --matrix ../../Trinity.gene.TMM.EXPR.matrix --samples ../../samples.txt -P 1e-3 -C 2

###TRINOTATE

/home/anabwire/bin/Trinotate-Trinotate-v4.0.0/Trinotate --create \
                          --db myTrinotate.sqlite \
                          --trinotate_data_dir /home/anabwire/data/aquality/TRINOTATE_DATA_DIR \
                          --use_diamond


##transdecode a two step process
####First run the TransDecoder step that identifies all long ORFs.

/home/anabwire/bin/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs -t ../Trinity.fasta

###Now, run the step that predicts which ORFs are likely to be coding.

/home/anabwire/bin/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict -t ../Trinity.fasta


nohup /home/anabwire/bin/Trinotate-Trinotate-v4.0.0/Trinotate --db myTrinotate.sqlite --init \
          --gene_trans_map ../Trinity.fasta.gene_trans_map \
          --transcript_fasta ../Trinity.fasta \
          --transdecoder_pep Trinity.fasta.transdecoder.pep >> trinotate1.out

nohup  /home/anabwire/bin/Trinotate-Trinotate-v4.0.0/Trinotate --db myTrinotate.sqlite --CPU 8 \
                     --transcript_fasta ../Trinity.fasta \
                     --transdecoder_pep Trinity.fasta.transdecoder.pep \
                     --trinotate_data_dir /home/anabwire/data/aquality/TRINOTATE_DATA_DIR \
                     --run "swissprot_blastp swissprot_blastx pfam EggnogMapper" \
                     --use_diamond >> trinotate2.out &

#generating trinotate report

Trinotate --db <sqlite.db> --report [ -E (default: 1e-5) ]
              [--pfam_cutoff DNC|DGC|DTC|SNC|SGC|STC (default: DNC=domain noise cutoff)]
              [--incl_pep]
              [--incl_trans]


/home/anabwire/bin/Trinotate-Trinotate-v4.0.0/Trinotate --db myTrinotate.sqlite --report > myTrinotate.tsv


## functional annotation

/home/anabwire/bin/Trinotate-Trinotate-v4.0.0/util/Trinotate_get_feature_name_encoding_attributes.pl \
                  myTrinotate.tsv  > annot_feature_map.txt

/home/anabwire/miniconda3/envs/trinity/bin/rename_matrix_feature_identifiers.pl \
                      ../Trinity.isoform.counts.matrix annot_feature_map.txt > Trinity_trans.counts.wAnnot.matrix

/home/anabwire/trinityrnaseq-v2.15.1/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl \
                      ../Trinity.isoform.counts.matrix annot_feature_map.txt > Trinity_trans.counts.wAnnot.matrix

###go ANALYSIS

nohup /home/anabwire/bin/Trinotate-Trinotate-v4.0.0/util/extract_GO_assignments_from_Trinotate_xls.pl \
                       --Trinotate_xls myTrinotate.tsv \
                       -G --include_ancestral_terms \
                       > go_annotations.txt >> go.out &
                       ##could no go ahead I need the file with the gene of interest so I propceeded to do the go analysis with the transcripts afterr DE analysis
/home/anabwire/trinityrnaseq-v2.15.1/util/misc/fasta_seq_length.pl  ../Trinity.fasta > Trinity.fasta.seq_lens

## create a gene_lengths file
/home/anabwire/trinityrnaseq-v2.15.1/util/misc/TPM_weighted_gene_length.py  \
         --gene_trans_map ../Trinity.fasta.gene_trans_map \
         --trans_lengths Trinity.fasta.seq_lens \
         --TPM_matrix ../Trinity.isoform.TMM.EXPR.matrix > Trinity.gene_lengths.txt

         /home/anabwire/trinityrnaseq-v2.15.1/util/misc/TPM_weighted_gene_length.py  \
                  --gene_trans_map ../Trinity.fasta.gene_trans_map \
                  --trans_lengths Trinity.fasta.seq_lens \
                  --TPM_matrix ../Trinity.gene.TMM.EXPR.matrix > Trinity.gene_lengths1.txt

##go enrichmentt ANALYSIS
/home/anabwire/trinityrnaseq-v2.15.1/Analysis/DifferentialExpression/run_GOseq.pl \
                      --factor_labeling  ../factor_labeling.txt \
                      --GO_assignments ../go_annotations.txt1 \
                      --lengths ../gene.lengths.txt \
                      --background /home/anabwire/data/aquality/DEgene/edgeR.1195828.dir/allgenes.txt

/home/anabwire/trinityrnaseq-v2.15.1/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl \
                        -R  ../edgeR.622143.dir/diffExpr.P1e-3_C2.matrix.RData --Ptree 60

/home/anabwire/bin/Trinotate-Trinotate-v4.0.0/util/Trinotate_get_feature_name_encoding_attributes.pl \
                                ../trinotate/myTrinotate.tsv > Trinotate_report.xls.name_mappings

####update your expression matrix to incorporate these new function-encoded feature identifiers
/home/anabwire/trinityrnaseq-v2.15.1/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl \
  ../Trinity.isoform.TMM.EXPR.matrix Trinotate_report.xls.name_mappings  > Trinity_trans.TMM.EXPR.annotated.matrix




## cleanup
## I removed the  corrected_reads from ./data to ./data/trinity/
