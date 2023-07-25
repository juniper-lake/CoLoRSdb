version 1.0

# Calculate VCF stats

import "../structs.wdl"

task small_variant_stats {
  input {
    File vcf
    Array[String] sample_ids

    File reference
    File non_diploid_regions

    RuntimeAttributes runtime_attributes
  }

  String vcf_basename = basename(vcf, ".gz")

  Int threads = 2
  Int reference_size = if (defined(reference)) then ceil(size(reference, "GB")) else 0
  Int disk_size = ceil((size(vcf, "GB") + reference_size) * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --help

    bcftools stats \
      --threads ~{threads - 1} \
      --samples ~{sep="," sample_ids} \
      --targets-file ^~{non_diploid_regions} \
      --fasta-ref ~{reference} \
      ~{vcf} \
    > ~{vcf_basename}.diploid_regions.stats.txt
  >>>

  output {
    File stats = "~{vcf_basename}.diploid_regions.stats.txt"
  }

  runtime {
    cpu: threads
    memory: "4 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:2f1bdbe671c7a5335ba75263347f81a48a0b56df523001aabdfe622fe825ed82"
  }
}

task structural_variant_stats {
  input {
    File vcf
    Array[String] sample_ids

    File non_diploid_regions

    RuntimeAttributes runtime_attributes
  }

  String vcf_basename = basename(vcf, ".gz")

  Int threads = 2
  Int disk_size = ceil((size(vcf, "GB")) * 1.5 + 20)

  command <<<
    set -euo pipefail

    bcftools --help

    for SAMPLE in ~{sep=" " sample_ids}; do
      bcftools query \
        --targets-file ^~{non_diploid_regions} \
        --exclude 'GT="miss"' \
        --samples $SAMPLE \
        --format '%SVTYPE\t[%GT]\n]' \
        ~{vcf} > ~{vcf_basename}.$SAMPLE.txt

      cat ~{vcf_basename}.$SAMPLE.txt \
        | datamash -s groupby 1,2 count 1 \
        | awk -v OFS="\t" -v sample=$SAMPLE '{print sample,$0}' \
        >> ~{vcf_basename}.diploid_regions.stats.tsv
    done
  >>>

  output {
    File stats = "~{vcf_basename}.diploid_regions.stats.tsv"
  }

  runtime {
    cpu: threads
    memory: "4 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:2f1bdbe671c7a5335ba75263347f81a48a0b56df523001aabdfe622fe825ed82"
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
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:2f1bdbe671c7a5335ba75263347f81a48a0b56df523001aabdfe622fe825ed82"
  }
}
