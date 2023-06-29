version 1.0

import "../../tasks/pbmm2.wdl" as Pbmm2
import "../../tasks/somalier.wdl" as Somalier
import "../../tasks/mosdepth.wdl" as Mosdepth

workflow cohort_align_qc {

  input {
    String cohort_id
    Array[Sample] samples

    Float max_sample_relatedness_qc
    Float min_movie_relatedness_qc

    ReferenceData reference

    RuntimeAttributes default_runtime_attributes
  }

  scatter (sample in samples) {
    # align each movie with pbmm2
    scatter (idx in range(length(sample.movies))) {    

      # get movie name from movie path but append "_movie1", "_movie2", etc in case file names
      # are not unique or don't follow HiFi naming convention
      String movie_name = sub(basename(sample.movies[idx]), "\\..*", "_movie~{idx}")
      call Pbmm2.pbmm2 {
          input: 
            movie = sample.movies[idx],
            out_prefix = "~{sample.sample_id}.~{movie_name}",
            sample_id = sample.sample_id,
            reference_name = reference.name,
            reference_fasta = reference.fasta.data,
            reference_index = reference.fasta.index,
            runtime_attributes = default_runtime_attributes
      }
      
      IndexData aligned_bam = {
        "data": pbmm2.aligned_bam,
        "index": pbmm2.aligned_bam_index
      }
      
      # extract somalier sites from aligned bam
      call Somalier.somalier_extract {
        input:
          sample_id = sample.sample_id,
          sample_prefix = "~{movie_name}_",
          bam = aligned_bam.data,
          bam_index = aligned_bam.index,
          reference_fasta = reference.fasta.data,
          reference_index = reference.fasta.index,
          somalier_sites_vcf = reference.somalier_sites_vcf,
          runtime_attributes = default_runtime_attributes
      }
    }

    # combine smrtcells stats for multiple movies
    call Pbmm2.combine_smrtcell_stats {
      input:
        sample_id = sample.sample_id,
        read_length_and_quality_tsvs = pbmm2.bam_stats,
        runtime_attributes = default_runtime_attributes
    }

    # check for sample swaps among movies
    call Somalier.somalier_relate_movies {
      input:
        sample_id = sample.sample_id,
        extracted_sites_per_movie = somalier_extract.extracted_sites_per_movie,
        runtime_attributes = default_runtime_attributes
    }

    Boolean pass_swap = if somalier_relate_movies.min_relatedness >= min_movie_relatedness_qc then true else false

    # merge aligned bams
    call Pbmm2.merge_bams {
      input:
        bams = pbmm2.aligned_bam,
        output_bam_name = "~{sample.sample_id}.~{reference.name}.bam",
        runtime_attributes = default_runtime_attributes
    }

    IndexData merged_bam = {
      "data": merge_bams.merged_bam, 
      "index": merge_bams.merged_bam_index
    }

    # calculate depth of merged bam
    call Mosdepth.mosdepth {
      input:
        aligned_bam = merged_bam.data,
        aligned_bam_index = merged_bam.index,
        runtime_attributes = default_runtime_attributes
    }
    
    String sample_id = sample.sample_id
    Int n_movies = length(sample.movies)
  }

  # somalier_relate_samples
  call Somalier.somalier_relate_samples {
    input:
      cohort_id = cohort_id,
      extracted_sites_per_sample = flatten(somalier_extract.extracted_sites_per_sample),
      sample_ids = sample_id,
      coverages = mosdepth.mean_coverage,
      max_pairwise_relatedness = max_sample_relatedness_qc,
      runtime_attributes = default_runtime_attributes
  }

  scatter (idx in range(length(samples))) {
    
    Boolean pass_sex = if somalier_relate_samples.inferred_sexes[idx] != "unknown" then true else false
    Boolean pass_relatedness = if somalier_relate_samples.qc_related_keep_drop[idx] == 'keep' then true else false
    Boolean qc_pass_combined = if pass_swap[idx] && pass_relatedness && pass_sex then true else false
    Float sample_relatedness_qc_threshold = max_sample_relatedness_qc
    Float movie_relatedness_qc_threshold = min_movie_relatedness_qc
    String min_movie_relatedness = "~{somalier_relate_movies.min_relatedness[idx]}"
    String n_relations = "~{somalier_relate_samples.n_relations[idx]}"
  }

  call summarize_qc {
    input:
      cohort_id = cohort_id,
      reference_name = reference.name,
      sample_ids = sample_id,
      movie_relatedness_qc_threshold = movie_relatedness_qc_threshold,
      sample_relatedness_qc_threshold = sample_relatedness_qc_threshold,
      n_movies = n_movies,
      qc_pass_sex = pass_sex,
      qc_pass_swap = pass_swap,
      qc_pass_relatedness = pass_relatedness,
      qc_pass_combined = qc_pass_combined,
      min_movie_relatedness = min_movie_relatedness,
      n_relations = n_relations,
      sex = somalier_relate_samples.inferred_sexes,
      coverage_mean = mosdepth.mean_coverage,
      read_count = combine_smrtcell_stats.read_count,
      unique_read_count = combine_smrtcell_stats.unique_read_count,
      read_quality_mean = combine_smrtcell_stats.read_quality_mean,
      read_quality_median = combine_smrtcell_stats.read_quality_median,
      read_quality_stdev = combine_smrtcell_stats.read_quality_stdev,
      read_length_mean = combine_smrtcell_stats.read_length_mean,
      read_length_median = combine_smrtcell_stats.read_length_median,
      read_length_stdev = combine_smrtcell_stats.read_length_stdev,
      runtime_attributes = default_runtime_attributes
  }

  output {
    Array[IndexData] aligned_bams = merged_bam
    Array[Boolean] qc_pass = qc_pass_combined
    Array[String] inferred_sexes = somalier_relate_samples.inferred_sexes
    File pairwise_relatedness = somalier_relate_samples.pairs
    File qc_summary_tsv = summarize_qc.quality_control_summary
  }
}


task summarize_qc {
  input {
    String cohort_id
    String reference_name
    Array[String] sample_ids
    Array[Float] movie_relatedness_qc_threshold
    Array[Float] sample_relatedness_qc_threshold
    Array[Boolean] qc_pass_swap
    Array[Boolean] qc_pass_sex
    Array[Boolean] qc_pass_relatedness
    Array[Boolean] qc_pass_combined
    Array[String] min_movie_relatedness
    Array[String] n_relations
    Array[String] sex
    Array[Int] n_movies
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
    
    paste <(echo -e "sample_id\n~{sep="\n" sample_ids}") \
      <(echo -e "movie_relatedness_qc_threshold\n~{sep="\n" movie_relatedness_qc_threshold}") \
      <(echo -e "sample_relatedness_qc_threshold\n~{sep="\n" sample_relatedness_qc_threshold}") \
      <(echo -e "qc_pass_swap\n~{sep="\n" qc_pass_swap}") \
      <(echo -e "qc_pass_relatedness\n~{sep="\n" qc_pass_relatedness}") \
      <(echo -e "qc_pass_sex\n~{sep="\n" qc_pass_sex}") \
      <(echo -e "qc_pass_combined\n~{sep="\n" qc_pass_combined}") \
      <(echo -e "min_movie_relatedness\n~{sep="\n" min_movie_relatedness}") \
      <(echo -e "n_relations\n~{sep="\n" n_relations}") \
      <(echo -e "sex\n~{sep="\n" sex}") \
      <(echo -e "n_movies\n~{sep="\n" n_movies}") \
      <(echo -e "coverage_mean\n~{sep="\n" coverage_mean}") \
      <(echo -e "read_count\n~{sep="\n" read_count}") \
      <(echo -e "unique_read_count\n~{sep="\n" unique_read_count}") \
      <(echo -e "read_quality_mean\n~{sep="\n" read_quality_mean}") \
      <(echo -e "read_quality_median\n~{sep="\n" read_quality_median}") \
      <(echo -e "read_quality_stdev\n~{sep="\n" read_quality_stdev}") \
      <(echo -e "read_length_mean\n~{sep="\n" read_length_mean}") \
      <(echo -e "read_length_median\n~{sep="\n" read_length_median}") \
      <(echo -e "read_length_stdev\n~{sep="\n" read_length_stdev}") \
      > ~{cohort_id}.quality_control_summary.tsv
  >>>

  output {
    File quality_control_summary = "~{cohort_id}.~{reference_name}.quality_control_summary.tsv"
  }
  runtime {
    cpu: 1
    memory: "1 GB"
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
