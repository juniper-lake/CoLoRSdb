version 1.0

# Merge pbsv SVs with jasminesv

import "../structs.wdl"

task jasminesv_merge_svs {
  input {
    Array[File] vcfs
    String output_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads = 16
  Int mem_gb = 64
  Int disk_size = ceil((size(vcfs, "GB")) * 3 + 20)

  command <<<
    set -euo pipefail

    # increase open file limit
    ulimit -Sn 65536

    jasmine --help

    jasmine \
      --output_genotypes \
      --allow_intrasample \
      --comma_filelist \
      threads=~{threads} \
      file_list=~{sep="," vcfs} \
      out_file=~{output_prefix}.vcf
  >>>

  output {
    File merged_vcf = "~{output_prefix}.vcf"
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/jasminesv@sha256:17e6ff136c32342e396f19a3bc83a0c84ae75e7018ade7d8c6976825e35d76ff"
  }
}
