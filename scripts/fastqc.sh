rule_rawfastqc:
  input:
    /home/anabwire/data/raw_data/{sample}.fastq.gz
  output:
    FASTQC/
  log:
          "logs/fastqc/{sample}.log"
  shell:
  "mkdir FASTQC",
  "fastqc -o {input}"

  rule fastqc:
      input:
          "reads/{sample}.fastq"
      output:
          html="qc/fastqc/{sample}.html",
          zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
      params: "--quiet"
      log:
          "logs/fastqc/{sample}.log"
      threads: 2
      wrapper:
          "/usr/local/bin/fastqc"
