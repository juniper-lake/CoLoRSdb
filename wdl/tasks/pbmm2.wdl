version 1.0

# Align HiFi reads fastq or bam with pbmm2

import "../structs.wdl"

task pbmm2 {

  input {
    File movie
    String sample_id
    
    String reference_name
    File reference_fasta
    File reference_index
    
    RuntimeAttributes runtime_attributes
    }

  String movie_name = sub(basename(movie), "\\..*", "")

  Int threads = 24
  Int mem_gb = ceil(threads * 4)
	Int disk_size = ceil((size(movie, "GB") + size(reference_fasta, "GB")) * 4 + 20)
  
  command {
    set -euo pipefail
    
    pbmm2 align \
      -num-threads ~{threads} \
      --sort-memory 4G \
      --sample ~{sample_id} \
      --log-level INFO \
      --preset CCS \
      --sort \
      --unmapped \
      ~{reference_fasta} \
      ~{movie} \
      ~{sample_id}.~{movie_name}.~{reference_name}.bam
  
    # movie stats
		extract_read_length_and_qual.py \
			~{movie} \
		> ~{sample_id}.~{movie_name}.read_length_and_quality.tsv
    }



  output {
    File aligned_bam = "~{sample_id}.~{movie_name}.~{reference_name}.bam"
    File aligned_bam_index = "~{sample_id}.~{movie_name}.~{reference_name}.bam.bai"
    File bam_stats = "~{sample_id}.~{movie_name}.read_length_and_quality.tsv"
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
