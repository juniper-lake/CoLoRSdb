version 1.0

struct ReferenceData {
  String name
  File fasta
  File fasta_index
  Array[String]+ chromosomes
  Array[String]+ autosomes
  File non_diploid_regions
  File tandem_repeat_bed
  File somalier_sites_vcf

  Array[File]? trgt_tandem_repeat_beds
  File? hificnv_exclude_bed
  File? hificnv_exclude_bed_index
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
