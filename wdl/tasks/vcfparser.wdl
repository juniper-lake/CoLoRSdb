version 1.0

# Postprocess VCFs by anonymizing, fixing ploidy, and zipping/indexing

import "../structs.wdl"

task postprocess_joint_vcf {
  input {
    File vcf
    String cohort_id
    Boolean anonymize_output
    Array[String]? sample_plus_sexes
    File? non_diploid_regions
    
    RuntimeAttributes runtime_attributes
  }

  String vcf_basename = basename(basename(vcf, ".gz"), ".vcf")
  String anonymize_prefix = if (anonymize_output) then "--anonymize_prefix ~{cohort_id}" else ""
  String outfile = if (anonymize_output) then "~{vcf_basename}.anonymized.vcf" else "~{vcf_basename}.ploidy_fixed.vcf"
  
  Int threads = 4
  Int disk_size = ceil((size(vcf, "GB")) * 3.5 + 20) 
  
  command <<<
    set -euo pipefail

    if ~{defined(sample_plus_sexes)}; then
      SAMPLE_SEXES="--sample_sexes ~{sep=' ' sample_plus_sexes}"
    else
      SAMPLE_SEXES=""
    fi
    
    postprocess_joint_vcf.py \
      ~{anonymize_prefix} \
      ~{"--non_diploid_regions " + non_diploid_regions} \
      $SAMPLE_SEXES \
      --outfile ~{outfile} \
      ~{vcf}

    bgzip --version

    bgzip \
      --threads ~{threads} \
      ~{outfile}

    tabix --version

    tabix \
      --preset vcf \
      ~{outfile}.gz   
  >>>

  output {
    File postprocessed_vcf = "~{outfile}.gz"
    File postprocessed_vcf_index = "~{outfile}.gz.tbi"
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
    docker: "~{runtime_attributes.container_registry}/vcfparser@sha256:f6e5b3c8fd0f0fdc8c67a00e084817d73d15971808f48aabb3e8a82b9c9c4f5f"
  }
}

task merge_trgt_vcfs {
  input {
    Array[File] trgt_vcfs
    String cohort_id
    String reference_name
    Boolean anonymize_output

    RuntimeAttributes runtime_attributes
  }

  String anonymize_prefix = if (anonymize_output) then "--anonymize_prefix " + cohort_id else ""
  String outfile = if (anonymize_output) then "~{cohort_id}.~{reference_name}.trgt.anonymized.vcf" else "~{cohort_id}.~{reference_name}.trgt.vcf"
  Int threads = 4
  Int disk_size = ceil((size(trgt_vcfs, "GB")) * 3.5 + 20)

  command {
    set -euo pipefail

    # increase open file limit
    ulimit -Sn 65536
    
    merge_trgt_vcfs.py \
      --outfile ~{outfile} \
      ~{anonymize_prefix} \
      ~{sep=" " trgt_vcfs} 

    bgzip \
      --threads ~{threads} \
      ~{outfile} \

    tabix \
      --preset vcf \
      ~{outfile}.gz
  }

  output {
    File merged_vcf = "~{outfile}.gz"
    File merged_vcf_index = "~{outfile}.gz.tbi"
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
    docker: "~{runtime_attributes.container_registry}/vcfparser@sha256:f6e5b3c8fd0f0fdc8c67a00e084817d73d15971808f48aabb3e8a82b9c9c4f5f"
  }
}
