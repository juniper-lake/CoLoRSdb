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
  Int threads = 4
  Int disk_size = ceil((size(bam, "GB") + size(reference_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    trgt \
      --karyotype ~{karyotype} \
      --genome ~{reference_fasta} \
      --repeats ~{tandem_repeat_bed} \
      --reads ~{bam} \
      --output-prefix ~{bam_basename}.trgt

    bcftools sort \
      --output-type z \
      --output ~{bam_basename}.trgt.sorted.vcf.gz \
      ~{bam_basename}.trgt.vcf.gz
  >>>

  output {
    File repeat_vcf = "~{bam_basename}.trgt.sorted.vcf.gz"
  }

  runtime {
    cpu: threads
    memory: "4 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/trgt@sha256:daeb2d091ccf412f9ba92229cda02f110b7ec7fa31f01453796170ed62664241"
  }
}
