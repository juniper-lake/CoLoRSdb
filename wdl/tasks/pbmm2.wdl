version 1.0

# Align HiFi reads fastq or bam with pbmm2

import "../structs.wdl"

task pbmm2 {
  input {
    File movie
    String out_prefix
    String sample_id

    String reference_name
    File reference_fasta
    File reference_index

    RuntimeAttributes runtime_attributes
    }

  Int threads = 24
  Int mem_gb = ceil(threads * 4)
  Int disk_size = ceil((size(movie, "GB") + size(reference_fasta, "GB")) * 4 + 20)

  command {
    set -euo pipefail

    pbmm2 align \
      --num-threads ~{threads} \
      --sort-memory 4G \
      --sample ~{sample_id} \
      --log-level INFO \
      --preset CCS \
      --sort \
      --unmapped \
      ~{reference_fasta} \
      ~{movie} \
      ~{out_prefix}.~{reference_name}.bam

    # movie stats
    extract_read_length_and_qual.py \
      ~{movie} \
    > ~{out_prefix}.read_length_and_quality.tsv
  }

  output {
    File aligned_bam = "~{out_prefix}.~{reference_name}.bam"
    File aligned_bam_index = "~{out_prefix}.~{reference_name}.bam.bai"
    File bam_stats = "~{out_prefix}.read_length_and_quality.tsv"
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
    docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:c9c9dcb98b80de8e877e3bca5c035622e15d5e68341f86bbcecfadb771389177"
  }
}

task combine_smrtcell_stats {
  input {
    String sample_id
    Array[File] read_length_and_quality_tsvs

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(read_length_and_quality_tsvs, "GB") + 20)

  command <<<
    set -euo pipefail

    cat ~{sep=' ' read_length_and_quality_tsvs} \
      | datamash count 1 countunique 1 median 2 mean 2 pstdev 2 median 3 mean 3 pstdev 3 \
      | awk '{OFS="\n"; $1=$1}1' \
      > ~{sample_id}.smrtcell_stats.txt

    sed '1q;d' ~{sample_id}.smrtcell_stats.txt > read_count.txt
    sed '2q;d' ~{sample_id}.smrtcell_stats.txt > read_count_unique.txt
    sed '3q;d' ~{sample_id}.smrtcell_stats.txt > read_length_median.txt
    sed '4q;d' ~{sample_id}.smrtcell_stats.txt > read_length_mean.txt
    sed '5q;d' ~{sample_id}.smrtcell_stats.txt > read_length_stdev.txt
    sed '6q;d' ~{sample_id}.smrtcell_stats.txt > read_quality_median.txt
    sed '7q;d' ~{sample_id}.smrtcell_stats.txt > read_quality_mean.txt
    sed '8q;d' ~{sample_id}.smrtcell_stats.txt > read_quality_stdev.txt
  >>>

  output {
    Int read_count = read_int("read_count.txt")
    Int unique_read_count = read_int("read_count_unique.txt")
    Int read_length_median = read_int("read_length_median.txt")
    Float read_length_mean = read_float("read_length_mean.txt")
    Float read_length_stdev = read_float("read_length_stdev.txt")
    Int read_quality_median = read_int("read_quality_median.txt")
    Float read_quality_mean = read_float("read_quality_mean.txt")
    Float read_quality_stdev = read_float("read_quality_stdev.txt")
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
    docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:c9c9dcb98b80de8e877e3bca5c035622e15d5e68341f86bbcecfadb771389177"
  }
}
