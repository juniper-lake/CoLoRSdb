version 1.0

# Quality control of BAMs using peek-a-bam.py

import "../structs.wdl"

task peek_a_bam {
  input {
    File bam

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(bam, "GB")) * 1.5 + 5)
  Int threads = 1

  command<<<
    set -euo pipefail

    # look at first N reads and get a summary of the BAM
    peek-a-bam.py ~{bam} \
      | jq -cr '. | [
      .file,
      (.samples | join(",")),
      .potentially_multiplexed,
      .demultiplexed,
      .ccs,
      .hifi,
      .kinetics,
      .base_modification,
      .aligned,
      .haplotagged
      ] | @tsv' > peekabam.tsv

    qc_pass_bam=true
    # if more than one sample, fail
    if [ $(cut -f2 peekabam_out.tsv | tr ',' '\n' | wc -l) -gt 1 ]; then
      qc_pass_bam=false;
    fi
    # if potentially multiplexed, fail
    if [ $(cut -f3 peekabam_out.tsv) = true ]; then
      qc_pass_bam=false;
    fi
    # if not demultiplexed, fail
    if [ $(cut -f4 peekabam_out.tsv) = false ]; then
      qc_pass_bam=false;
    fi
    # if not ccs, fail
    if [ $(cut -f5 peekabam_out.tsv) = false ]; then
      qc_pass_bam=false;
    fi
    # if not hifi, fail
    if [ $(cut -f6 peekabam_out.tsv) = false ]; then
      qc_pass_bam=false;
    fi

    echo $qc_pass_bam > qc_pass_bam.txt
  >>>

  output {
    Array[String] bam_qc_tsv = read_tsv("peekabam.tsv")[0]
    Boolean qc_pass_bam = read_boolean("qc_pass_bam.txt")
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
    docker: "~{runtime_attributes.container_registry}/peekabam@sha256:380820974acb6bf29d624f43c42ae8d2889e38cf8d943c14b37cc176031eab1c"
  }
}
