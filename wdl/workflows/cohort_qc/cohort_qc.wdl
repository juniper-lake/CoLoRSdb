version 1.0

import "../../tasks/somalier.wdl" as Somalier

workflow cohort_qc {

  input {
    String cohort_id

    Float max_pairwise_relatedness

    Array[String] sample_ids
    Array[File] extracted_somalier_sites
    Array[Int] n_movies
    Array[Boolean] qc_pass_sex
    Array[Boolean] qc_pass_swap
    Array[String] sex
    Array[Float] min_movie_relatedness
    Array[Float] coverage
    Array[Int] read_count
    Array[Int] unique_read_count
    Array[Float] read_quality_mean
    Array[Float] read_quality_median
    Array[Float] read_quality_stdev
    Array[Float] read_length_mean
    Array[Float] read_length_median
    Array[Float] read_length_stdev

    RuntimeAttributes default_runtime_attributes
  }

  # somalier_relate_samples
  call Somalier.somalier_relate_samples {
    input:
      cohort_id = cohort_id,
      extracted_somalier_sites = extracted_somalier_sites,
      sample_ids = sample_ids,
      coverages = coverage,
      max_pairwise_relatedness = max_pairwise_relatedness,
      runtime_attributes = default_runtime_attributes
  }

  scatter (keep_drop in somalier_relate_samples.qc_related_keep_drop[1]) {
    Boolean pass_relatedness = if keep_drop == 'keep' then true else false
  }

  scatter (idx in range(length(sample_ids))){
    Boolean qc_pass_combined = if qc_pass_sex[idx] && qc_pass_swap[idx] && pass_relatedness[idx] then true else false
  }

  call summarize_qc {
    input:
      cohort_id = cohort_id,
      sample_ids = sample_ids,
      n_movies = n_movies,
      qc_pass_sex = qc_pass_sex,
      qc_pass_swap = qc_pass_swap,
      qc_pass_relatedness = pass_relatedness,
      qc_pass_combined = qc_pass_combined,
      sex = sex,
      min_movie_relatedness = min_movie_relatedness,
      coverage_mean = coverage,
      read_count = read_count,
      unique_read_count = unique_read_count,
      read_quality_mean = read_quality_mean,
      read_quality_median = read_quality_median,
      read_quality_stdev = read_quality_stdev,
      read_length_mean = read_length_mean,
      read_length_median = read_length_median,
      read_length_stdev = read_length_stdev,
      runtime_attributes = default_runtime_attributes
  }

  output {
    Array[Boolean] qc_pass = qc_pass_combined
    File pairwise_relatedness = somalier_relate_samples.pairs
    File qc_summary = summarize_qc.quality_control_summary
  }
}


task summarize_qc {
  input {
    String cohort_id
    Array[String] sample_ids
    Array[Int] n_movies
    Array[Boolean] qc_pass_sex
    Array[Boolean] qc_pass_swap
    Array[Boolean] qc_pass_relatedness
    Array[Boolean] qc_pass_combined
    Array[String] sex
    Array[Float] min_movie_relatedness
    Array[Float] coverage_mean
    Array[Float] read_count
    Array[Float] unique_read_count
    Array[Float] read_quality_mean
    Array[Float] read_quality_median
    Array[Float] read_quality_stdev
    Array[Float] read_length_mean
    Array[Float] read_length_median
    Array[Float] read_length_stdev

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 20

  command <<<
    set -euo pipefail
    echo -e "sample_id\tn_movies\tqc_pass_sex\tqc_pass_swap\tqc_pass_relatedness\tqc_pass_combined\tsex\tmin_movie_relatedness\tcoverage_mean\tread_count\tunique_read_count\tread_quality_mean\tread_quality_median\tread_quality_stdev\tread_length_mean\tread_length_median\tread_length_stdev" \
      > ~{cohort_id}.quality_control_summary.tsv
    
    paste <(echo -e ~{sep="\n" sample_ids}) \
      <(echo -e ~{sep="\n" n_movies}) \
      <(echo -e ~{sep="\n" qc_pass_sex}) \
      <(echo -e ~{sep="\n" qc_pass_swap}) \
      <(echo -e ~{sep="\n" qc_pass_relatedness}) \
      <(echo -e ~{sep="\n" qc_pass_combined}) \
      <(echo -e ~{sep="\n" sex}) \
      <(echo -e ~{sep="\n" min_movie_relatedness}) \
      <(echo -e ~{sep="\n" coverage_mean}) \
      <(echo -e ~{sep="\n" read_count}) \
      <(echo -e ~{sep="\n" unique_read_count}) \
      <(echo -e ~{sep="\n" read_quality_mean}) \
      <(echo -e ~{sep="\n" read_quality_median}) \
      <(echo -e ~{sep="\n" read_quality_stdev}) \
      <(echo -e ~{sep="\n" read_length_mean}) \
      <(echo -e ~{sep="\n" read_length_median}) \
      <(echo -e ~{sep="\n" read_length_stdev}) \
      >> ~{cohort_id}.quality_control_summary.tsv
  >>>

  output {
    File quality_control_summary = "~{cohort_id}.quality_control_summary.tsv"
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
    docker: "ubuntu:20.04"
  }
}