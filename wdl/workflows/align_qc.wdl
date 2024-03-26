version 1.0

# Align and perform QC on each sample in cohort

import "backend_configuration.wdl" as BackendConfiguration
import "../tasks/pbmm2.wdl" as Pbmm2
import "../tasks/somalier.wdl" as Somalier
import "../tasks/mosdepth.wdl" as Mosdepth
import "../tasks/samtools.wdl" as Samtools
import "../tasks/utils.wdl" as Utils

workflow align_qc {
  input {
    String cohort_id
    Array[String] sample_ids
    Array[Array[File]] sample_movies

    ReferenceData reference

    RuntimeAttributes default_runtime_attributes
  }

  # these inputs are not in the input section because should not be changed for data contributing to CoLoRSdb
  Float max_sample_relatedness_qc = 0.125
  Float min_movie_relatedness_qc = 0.7

  scatter (sample_idx in range(length(sample_ids))) {
    String sample_id = sample_ids[sample_idx]
    Array[File] movies = sample_movies[sample_idx]

    scatter (movie_idx in range(length(movies))) {
      # ignore FASTQs
      if (basename(movies[movie_idx], ".bam") != basename(movies[movie_idx])) {
        File bam = movies[movie_idx]
      }
      if (!(basename(movies[movie_idx], ".bam") != basename(movies[movie_idx]))) {
        File fastq = movies[movie_idx]
      }
    }

    if (length(select_all(bam)) > 0) {
      # remove alignment info and unnecessary tags from bam s
      call Samtools.reset_bams_and_qc {
        input:
          sample_id = sample_id,
          bams = select_all(bam),
          runtime_attributes = default_runtime_attributes
      }
    }

    # check if bams pass QC, if no bams are present, set to true
    Boolean qc_pass_bam = select_first([reset_bams_and_qc.qc_pass, true])

    Array[File] movies_with_bams_reset = flatten([select_first([reset_bams_and_qc.reset_bams,[]]), select_all(fastq)])

    # align each movie with pbmm2
    scatter (movie_idx in range(length(movies_with_bams_reset))) {
      # get movie name from movie path but append "_movie1", "_movie2", etc in case file names
      # are not unique or don't follow HiFi naming convention
      String movie_name = sub(basename(movies_with_bams_reset[movie_idx]), "\\..*", "_movie~{movie_idx}")

      # align each movie
      call Pbmm2.pbmm2 {
          input:
            movie = movies_with_bams_reset[movie_idx],
            out_prefix = "~{sample_id}.~{movie_name}",
            sample_id = sample_id,
            reference_name = reference.name,
            reference_fasta = reference.fasta,
            reference_index = reference.fasta_index,
            runtime_attributes = default_runtime_attributes
      }
    }

    # combine smrtcells stats for multiple movies
    call Pbmm2.combine_smrtcell_stats {
      input:
        sample_id = sample_id,
        read_length_and_quality_tsvs = pbmm2.bam_stats,
        runtime_attributes = default_runtime_attributes
    }

    # check for sample swaps among movies
    call Somalier.somalier_sample_swap {
      input:
        sample_id = sample_id,
        bams = pbmm2.aligned_bam,
        bam_indexes = pbmm2.aligned_bam_index,
        movie_names = movie_name,
        reference_fasta = reference.fasta,
        reference_index = reference.fasta_index,
        somalier_sites_vcf = reference.somalier_sites_vcf,
        runtime_attributes = default_runtime_attributes
    }

    Boolean pass_swap = if somalier_sample_swap.min_relatedness >= min_movie_relatedness_qc then true else false

    # merge aligned bams
    call Samtools.merge_bams {
      input:
        bams = pbmm2.aligned_bam,
        output_bam_name = "~{sample_id}.~{reference.name}.bam",
        runtime_attributes = default_runtime_attributes
    }

    # extract somalier sites from merged bam
    call Somalier.somalier_extract as somalier_extract_merged {
      input:
        sample_id = sample_id,
        bam = merge_bams.merged_bam,
        bam_index = merge_bams.merged_bam_index,
        reference_fasta = reference.fasta,
        reference_index = reference.fasta_index,
        somalier_sites_vcf = reference.somalier_sites_vcf,
        runtime_attributes = default_runtime_attributes
    }

    # calculate depth of merged bam
    call Mosdepth.mosdepth {
      input:
        aligned_bam = merge_bams.merged_bam,
        aligned_bam_index = merge_bams.merged_bam_index,
        runtime_attributes = default_runtime_attributes
    }

    Int n_movies = length(movies)
  }

  # somalier_relate_samples
  call Somalier.somalier_relate_samples {
    input:
      cohort_id = cohort_id,
      reference_name = reference.name,
      extracted_sites = somalier_extract_merged.extracted_sites,
      sample_ids = sample_ids,
      coverages = mosdepth.mean_coverage,
      max_pairwise_relatedness = max_sample_relatedness_qc,
      runtime_attributes = default_runtime_attributes
  }

  scatter (idx in range(length(sample_ids))) {
    Boolean pass_sex = if somalier_relate_samples.inferred_sexes[idx] != "unknown" then true else false
    Boolean pass_relatedness = if somalier_relate_samples.qc_related_keep_drop[idx] == 'keep' then true else false
    Boolean qc_pass_combined = if pass_swap[idx] && pass_relatedness && pass_sex && qc_pass_bam[idx] then true else false
    Float sample_relatedness_qc_threshold = max_sample_relatedness_qc
    Float movie_relatedness_qc_threshold = min_movie_relatedness_qc

  }

  call Utils.summarize_qc {
    input:
      cohort_id = cohort_id,
      reference_name = reference.name,
      sample_ids = sample_ids,
      movie_relatedness_qc_threshold = movie_relatedness_qc_threshold,
      sample_relatedness_qc_threshold = sample_relatedness_qc_threshold,
      n_movies = n_movies,
      qc_pass_sex = pass_sex,
      qc_pass_swap = pass_swap,
      qc_pass_relatedness = pass_relatedness,
      qc_pass_bams = qc_pass_bam,
      qc_pass_combined = qc_pass_combined,
      min_movie_relatedness = somalier_sample_swap.min_relatedness,
      min_movie_relatedness_n_sites = somalier_sample_swap.min_relatedness_n_sites,
      n_relations = somalier_relate_samples.n_relations,
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
      peek_a_bam_tsv = flatten(select_all(reset_bams_and_qc.peekabam_tsv)),
      runtime_attributes = default_runtime_attributes
  }

  output {
    Array[String] sexes = somalier_relate_samples.inferred_sexes
    Array[Boolean] qc_pass = qc_pass_combined
    Array[File] aligned_bams = merge_bams.merged_bam
    Array[File] aligned_bam_indexes = merge_bams.merged_bam_index

    File somalier_pairs = somalier_relate_samples.pairs
    File somalier_samples = somalier_relate_samples.samples
    File sample_qc_summary_tsv = summarize_qc.sample_qc_tsv
    File bam_qc_summary_tsv = summarize_qc.bam_qc_tsv
  }
}
