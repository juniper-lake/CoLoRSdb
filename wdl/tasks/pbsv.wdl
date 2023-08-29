version 1.0

# Call structural variants with pbsv

import "../structs.wdl"

task pbsv_discover {
  input {
    File bam
    File bam_index

    String region
    File tandem_repeat_bed

    RuntimeAttributes runtime_attributes
    }

  String prefix = basename(bam, ".bam")

  Int disk_size = ceil((size(bam, "GB") + size(tandem_repeat_bed, "GB")) * 2.5 + 20)

  command<<<
    set -euo pipefail

    pbsv discover \
      --hifi \
      --region ~{region} \
      --tandem-repeats ~{tandem_repeat_bed} \
      ~{bam} \
      ~{prefix}.~{region}.svsig.gz
    >>>

  output {
    File svsig = "~{prefix}.~{region}.svsig.gz"
  }

  runtime {
    cpu: 2
    memory: "8 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/pbsv@sha256:616c2f7fa6832ea5493b4b1f0eaa62f286be9ec04d97d888d02e762a73b1737c"
  }
}

task pbsv_call {
  input {
    String sample_id
    Array[File] svsigs
    Int sample_count
    String region

    String reference_name
    File reference_fasta
    File reference_index

    Int? mem_gb

    RuntimeAttributes runtime_attributes
  }
  
  Int default_mem_gb = if sample_count > 20 then 128 else 64
  Int runtime_mem_gb = select_first([mem_gb, default_mem_gb])
  Int threads = 8
  Int disk_size = ceil((size(svsigs, "GB") + size(reference_fasta, "GB")) * 2 + 20)

  command<<<
    set -euo pipefail

    # increase open file limit
    ulimit -Sn 65536
    
    pbsv call \
      --log-level INFO \
      --hifi \
      --min-sv-length 20 \
      --types DEL,INS,INV \
      --num-threads ~{threads} \
      ~{reference_fasta} \
      ~{sep=" " svsigs} \
      ~{sample_id}.~{reference_name}.~{region}.pbsv.vcf
  >>>

  output {
    File vcf = "~{sample_id}.~{reference_name}.~{region}.pbsv.vcf"
    File log = stdout()
  }
  
  runtime {
    cpu: threads
    memory: "~{runtime_mem_gb} GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/pbsv@sha256:616c2f7fa6832ea5493b4b1f0eaa62f286be9ec04d97d888d02e762a73b1737c"
  }
}
