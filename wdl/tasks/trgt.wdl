version 1.0

# Genotype tandem repeats with TRGT

import "../structs.wdl"

task trgt {
  input {
    String sex

    File bam
    File bam_index

    File reference_fasta
    File reference_index
    File tandem_repeat_bed

    RuntimeAttributes runtime_attributes
  }

  String karyotype = if sex == "male" then "XY" else "XX"
  String bam_basename = basename(bam, ".bam")
  String bed_basename = basename(tandem_repeat_bed, ".bed")
  Int threads = 16
  Int mem_gb = 24
  Int disk_size = ceil((size(bam, "GB") + size(reference_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    trgt \
      --threads ~{threads} \
      --karyotype ~{karyotype} \
      --genome ~{reference_fasta} \
      --repeats ~{tandem_repeat_bed} \
      --reads ~{bam} \
      --output-prefix ~{bam_basename}.~{bed_basename}

    bcftools sort \
      --max-mem ~{mem_gb}G \
      --output-type z \
      --output ~{bam_basename}.~{bed_basename}.sorted.vcf.gz \
      ~{bam_basename}.~{bed_basename}.vcf.gz
  >>>

  output {
    File repeat_vcf = "~{bam_basename}.~{bed_basename}.sorted.vcf.gz"
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/trgt@sha256:7c759bdcbc07f4c5d812989fb0e06c3479c63622e7efc2d7299804e039f33610"
  }
}
