version 1.0

import "../structs.wdl"

task sniffles_discover {

  input {
    String sample_id
    File bam
    File bam_index

    File reference_fasta
    File reference_index
    File tandem_repeat_bed

    RuntimeAttributes runtime_attributes
  }

  String bam_basename = basename(bam, ".bam")
  String output_filename = "~{bam_basename}.snf"

  Int threads = 8
  Int mem_gb = 4 * threads
  Int disk_size = ceil(2.5 * (size(bam, "GB") + size(reference_fasta, "GB"))) + 20

  command {
    set -euo pipefail
    
    sniffles \
      --threads ~{threads} \
      --sample-id ~{sample_id} \
      --reference ~{reference_fasta} ~{"--tandem-repeats " + tandem_repeat_bed} \
      --minsvlen 20 \
      --snf ~{output_filename} \
      --input ~{bam}
  }

  output {
    File snf = output_filename
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
    docker: "~{runtime_attributes.container_registry}/sniffles:2.0.7"
  }
}


task sniffles_call {

  input {
    String sample_id
    Array[File] snfs

    String reference_name
    File reference_fasta
    File reference_index
    File tr_bed

    RuntimeAttributes runtime_attributes
  }

  String output_filename = "~{sample_id}.~{reference_name}.sniffles.vcf"

  Int threads = 8
  Int mem_gb = 4 * threads
  Int disk_size = ceil(2.5 * (size(snfs, "GB") + size(reference_fasta, "GB"))) + 20

  command {
    set -euo pipefail
    
    sniffles \
      --threads ~{threads} \
      --reference ~{reference_fasta} ~{"--tandem-repeats " + tr_bed} \
      --minsvlen 20 \
      --vcf ~{output_filename} \
      --input ~{sep=' ' snfs}
  }

  output {
    File vcf = output_filename
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
    docker: "~{runtime_attributes.container_registry}/sniffles:2.0.7"
  }
}
