version 1.0

struct IndexData {
  File data
  File index
}

struct Sample {
  String sample_id
  Array[File] movies
}

struct ReferenceData {
  String name
  IndexData fasta
  Array[String] chromosomes
  File non_diploid_regions
  File tandem_repeat_bed
  File somalier_sites_vcf

  File? trgt_tandem_repeat_bed
  IndexData? hificnv_exclude_bed
  File? hificnv_expected_bed_male
  File? hificnv_expected_bed_female
  File? peddy_sites
  File? peddy_bin
}

struct RuntimeAttributes {
  # The number of times to retry a task that fails due to preemption
  Int preemptible_tries
  # The number of times to retry a task that fails due a to nonzero return code
  Int max_retries

  String zones
  String queue_arn
  String container_registry
}
