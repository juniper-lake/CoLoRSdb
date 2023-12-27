version 1.0

# Calculate summary stats using mosdepth

import "../structs.wdl"

task mosdepth {
  input {
    File aligned_bam
    File aligned_bam_index

    RuntimeAttributes runtime_attributes
  }

  String prefix = basename(aligned_bam, ".bam")
  Int threads = 4
  Int disk_size = ceil(size(aligned_bam, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    mosdepth \
      --threads ~{threads - 1} \
      --by 500 \
      --no-per-base \
      --use-median \
      ~{prefix} \
      ~{aligned_bam}

    awk '$1=="total" { print $4 }' ~{prefix}.mosdepth.summary.txt > ~{prefix}.mean_coverage.txt
  >>>

  output {
    File summary = "~{prefix}.mosdepth.summary.txt"
    File region_bed = "~{prefix}.regions.bed.gz"
    Float mean_coverage = read_float("~{prefix}.mean_coverage.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/mosdepth@sha256:dc0331405d114151e65c4e1c4185b5e290141f7e2acebc2f18a9adefe908cfea"
    cpu: threads
    memory: "4 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " LOCAL"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
  }
}
