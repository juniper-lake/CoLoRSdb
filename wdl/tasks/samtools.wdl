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
    docker: "~{runtime_attributes.container_registry}/samtools@sha256:3e2f9f7761c9f704a2c5c0226f061f1a2f186e68c9001f05f0e424321f44aaa7"
  }
}

task reset_bams_and_qc {
  input {
    String sample_id
    Array[File]+ bams

    RuntimeAttributes runtime_attributes
  }

  Int threads = 8
  Int disk_size = ceil(size(bams, "GB") * 2.5 + 20)

  command <<<
    set -euo pipefail

    samtools --version

    qc_pass_bam=true

    for BAM in ~{sep=" " bams}; do
      # if not a fastq
      if [[ $BAM == *.bam ]]; then

        BASE=$(basename $BAM .bam)

        samtools reset \
          --remove-tag mc,mg,mi,rm,HP,PS,fi,ri,fp,rp,fn,rn \
          --threads ~{threads - 1} \
          -o ${BASE}.reset.bam \
          $BAM

        # look at first N reads and get a summary of the BAM
        peek-a-bam.py $BAM \
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
          ] | @tsv' > peekabam_out.tsv

        # if more than one sample, fail
        if [ $(cut -f2 peekabam_out.tsv | tr ',' '\n' | wc -l) -gt 1 ]; then
          qc_pass_bam=false;
        fi
        # if potentially multiplexed, fail
        if [ $(cut -f3 peekabam_out.tsv) = true ]; then
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
        awk -v OFS='\t' '{print "~{sample_id}", $0}' peekabam_out.tsv >> peekabam_with_ids.tsv
      fi
    done

    echo $qc_pass_bam > qc_pass_bam.txt

  >>>

  output {
    Array[File] reset_bams = glob("*.reset.bam")
    Array[Array[String]] peekabam_tsv = read_tsv("peekabam_with_ids.tsv")
    Boolean qc_pass = read_boolean("qc_pass_bam.txt")
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
    docker: "~{runtime_attributes.container_registry}/samtools@sha256:3e2f9f7761c9f704a2c5c0226f061f1a2f186e68c9001f05f0e424321f44aaa7"
  }
}
