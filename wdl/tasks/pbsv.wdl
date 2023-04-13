version 1.0

import "../structs.wdl"

task pbsv_discover {

  input {
    File bam
    File bam_index
    File tr_bed

    RuntimeAttributes runtime_attributes
    }

  String prefix = basename(bam, ".bam")
  Int threads = 16
  Int mem_gb = 4 * threads
  Int disk_size = ceil((size(bam, "GB") + size(tr_bed, "GB")) * 2.5 + 20)

  command<<<
    set -euo pipefail

    samtools view -H ~{bam} \
    | grep '^@SQ' \
    | cut -f2 \
    | cut -d':' -f2 \
    | parallel --jobs ~{threads} \
      pbsv discover \
        --log-level INFO \
        --hifi \
        --region {} \
        --tandem-repeats ~{tr_bed} \
        ~{bam} \
        ~{prefix}.{}.svsig.gz
    >>>

  output {
    Array[File] svsigs = glob("~{prefix}.*.svsig.gz")
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
    docker: "~{runtime_attributes.container_registry}/pbsv:2.9.0"
  }
}


task pbsv_call {

  input {
    String sample_id
    Array[File] svsigs

    String reference_name
    File reference_fasta
    File reference_index

    RuntimeAttributes runtime_attributes

  }
  
  Int threads = 8
  String output_filename = "~{sample_id}.~{reference_name}.pbsv.vcf"
  Int disk_size = ceil((size(svsigs, "GB") + size(reference_fasta, "GB")) * 2 + 20)

  command<<<
    set -euo pipefail

    pbsv call \
      --hifi \
      --min-sv-length 20 \
      --num-threads ~{threads} \
      ~{reference_fasta} \
      ~{sep=" " svsigs} \
      ~{output_filename}
  >>>

  output {
    File vcf = output_filename
  }
  
  runtime {
    cpu: threads
    memory: "64 GB"
		disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/pbsv:2.9.0"
  }
}
