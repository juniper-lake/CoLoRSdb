version 1.0

import "../structs.wdl"

task pbmm2_align {

  input {
    File movie
    String sample_id
    
    String reference_name
    File reference_fasta
    File reference_index
    
    RuntimeAttributes runtime_attributes
    }

  String movie_name = sub(basename(movie), "\\..*", "")
  String output_bam = "~{sample_id}.~{movie_name}.~{reference_name}.bam"

  Int threads = 24
  Int mem_gb = ceil(threads * 1.5)
	Int disk_size = ceil((size(movie, "GB") + size(reference_fasta, "GB")) * 4 + 20)
  
  command {
    set -o pipefail
    pbmm2 align \
      --sample ~{sample_id} \
      --log-level INFO \
      --preset CCS \
      --sort \
      --unmapped \
      -c 0 -y 70 \
      -j ~{threads} \
      ~{reference_fasta} \
      ~{movie} \
      ~{output_bam}
    }

  output {
    File bam = output_bam
    File bam_index = "~{output_bam}.bai"
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
    docker: "~{runtime_attributes.container_registry}/pbmm2:1.10.0"
  }
}
