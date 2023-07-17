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
      ~{prefix}.{region}.svsig.gz
    >>>

  output {
    File svsig = "~{prefix}.{region}.svsig.gz"
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
    docker: "~{runtime_attributes.container_registry}/pbsv:2.9.0"
  }
}

task pbsv_call {
  input {
    String sample_id
    Array[File] svsigs
    Int? sample_count
    String region

    String reference_name
    File reference_fasta
    File reference_index

    Int? pbsv_call_mem_gb

    RuntimeAttributes runtime_attributes
  }
  
  Int threads = 8
  Int default_mem_gb = if select_first([sample_count, 1]) > 3 then 96 else 64
  Int mem_gb = select_first([pbsv_call_mem_gb, default_mem_gb])
  Int disk_size = ceil((size(svsigs, "GB") + size(reference_fasta, "GB")) * 2 + 20)

  command<<<
    set -euo pipefail

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

task concat_vcfs {
  input {
    Array[File] vcfs
    String output_vcf_name

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs[0], "GB") * length(vcfs) * 2 + 20)

  command {
    set -euo pipefail

    if [[ "~{length(vcfs)}" -eq 1 ]]; then
      mv ~{vcfs[0]} ~{output_vcf_name}
    else
      # will throw error if vcfs are not in correct order
      # could --allow-overlaps but would require sort after and that's compute intensive
      bcftools concat \
        -o ~{output_vcf_name} \
        ~{sep=' ' vcfs}
    fi
  }

  output {
    File concatenated_vcf = output_vcf_name
  }

  runtime {
    cpu: 1
    memory: "2 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/bcftools:1.14"
  }
}
