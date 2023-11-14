version 1.0

# Utilities for samtools

import "../structs.wdl"

task merge_bams {
  input {
    Array[File] bams
    String output_bam_name

    RuntimeAttributes runtime_attributes
  }

  Int threads = 8
  Int disk_size = ceil(size(bams, "GB") * 2 + 20)

  command {
    set -euo pipefail

    samtools --version

    if [[ "~{length(bams)}" -eq 1 ]]; then
      cp ~{bams[0]} ~{output_bam_name}
    else
      samtools merge \
        -@ ~{threads - 1} \
        -o ~{output_bam_name} \
        ~{sep=' ' bams}
    fi

    samtools index ~{output_bam_name}
  }

  output {
    File merged_bam = output_bam_name
    File merged_bam_index = "~{output_bam_name}.bai"
  }

  runtime {
    cpu: threads
    memory: "1 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/samtools@sha256:a21de893caa867dbbedf5aa9403a6b5c1fc04c64bb0f8523a874c93f1b0edf08"
  }
}

task reset_aligned_bam {
  input {
    File bam

    RuntimeAttributes runtime_attributes
  }

  String bam_basename = basename(bam, ".bam")

  Int threads = 8
  Int disk_size = ceil(size(bam, "GB") * 2.5 + 20)

  command {
    set -euo pipefail

    samtools --version

    samtools reset \
      --remove-tag mc,mg,mi,rm,HP,fi,ri,fp,rp,fn,rn \
      --threads ~{threads - 1} \
      -o ~{bam_basename}.reset.bam \
      ~{bam}
  }

  output {
    File? unaligned_bam = "~{bam_basename}.reset.bam"
  }

  runtime {
    cpu: threads
    memory: "1 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/samtools@sha256:a21de893caa867dbbedf5aa9403a6b5c1fc04c64bb0f8523a874c93f1b0edf08"
  }
}
