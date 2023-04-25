version 1.0

import "../structs.wdl"

task hiphase {
  
  input {
    String cohort_id

    Array[String] sample_ids
    Array[File] aligned_bams
    Array[File] aligned_bam_indexes

    File deepvariant_vcf
    File? deepvariant_vcf_index
    File pbsv_vcf
    File? pbsv_vcf_index

    File reference_fasta
    File reference_index

	  RuntimeAttributes runtime_attributes
  }

  Int threads = 16
  Int mem_gb = 4 * threads
  Int disk_size = ceil(2.5 * (size(aligned_bams, "GB") + size(reference_fasta, "GB") + size(deepvariant_vcf, "GB") + size(pbsv_vcf, "GB"))) + 20


  String deepvariant_output_vcf_name = "~{basename(basename(deepvariant_vcf, '.gz'), '.vcf')}.hiphase.vcf.gz"
  String pbsv_output_vcf_name = "~{basename(basename(pbsv_vcf, '.gz'), '.vcf')}.hiphase.vcf.gz"

  command <<<
    set -euo pipefail

    hiphase \
      --threads ~{threads} \
      --reference ~{reference_fasta} \
      --global-realignment-cputime 300 \
      --bam ~{sep=" --bam " aligned_bams} \
      --sample-name ~{sep=" --sample-name " sample_ids} \
      --vcf ~{deepvariant_vcf} \
      --output-vcf ~{deepvariant_output_vcf_name} \
      --vcf ~{pbsv_vcf} \
      --output-vcf ~{pbsv_output_vcf_name} \
      --stats-file ~{cohort_id}.hiphase.stats.tsv \
      --blocks-file ~{cohort_id}.hiphase.blocks.tsv \
      --summary-file ~{cohort_id}.hiphase.summary.tsv
  >>>

  output {
    File deepvariant_output_vcf = deepvariant_output_vcf_name
    File pbsv_output_vcf = pbsv_output_vcf_name
    File stats = "~{cohort_id}.hiphase.stats.tsv"
    File blocks = "~{cohort_id}.hiphase.blocks.tsv"
    File summary = "~{cohort_id}.hiphase.summary.tsv"
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
    docker: "~{runtime_attributes.container_registry}/hiphase:0.8.0"
  }
}