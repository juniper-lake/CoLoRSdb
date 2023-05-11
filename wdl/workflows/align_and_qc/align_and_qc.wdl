version 1.0

import "../../tasks/pbmm2.wdl" as Pbmm2
import "../../tasks/somalier.wdl" as Somalier
import "../../tasks/utils.wdl" as Utils
import "../../tasks/mosdepth.wdl" as Mosdepth


workflow align_and_qc {

  input {
    String cohort_id
    Array[Sample] samples

    ReferenceData reference
    Float min_relatedness_sample_swap

    RuntimeAttributes default_runtime_attributes
  }

  scatter (sample in samples) {

    String sample_id = sample.sample_id
    Int n_movies = length(sample.movies)

    # align each movie with pbmm2
    scatter (idx in range(length(sample.movies))) {     
      call Pbmm2.pbmm2 {
          input: 
            movie = sample.movies[idx],
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

      # get movie name from movie path
      String movie_name = sub(basename(sample.movies[idx]), "\\..*", "")

      # extract somalier sites from aligned bam
      call Somalier.somalier_extract {
        input:
          sample_id = sample.sample_id,
          sample_prefix = movie_name + "_",
          bam = aligned_bam.data,
          bam_index = aligned_bam.index,
          reference_fasta = reference.fasta.data,
          reference_index = reference.fasta.index,
          somalier_sites_vcf = reference.somalier_sites_vcf,
          runtime_attributes = default_runtime_attributes
      }
    }

    # combine smrtcells stats
    call Utils.combine_smrtcell_stats {
      input:
        sample_id = sample.sample_id,
        read_length_and_quality_tsvs = pbmm2.bam_stats,
        runtime_attributes = default_runtime_attributes
    }

    # check for sample swaps among movies and confirm gender of sample
    call Somalier.somalier_relate_movies {
      input:
        sample_id = sample.sample_id,
        extracted_somalier_sites = somalier_extract.extracted_sites_per_movie,
        runtime_attributes = default_runtime_attributes
    }

    Boolean pass_sex = if somalier_relate_movies.inferred_sex != "UNKNOWN" then true else false
    Boolean pass_swap = if somalier_relate_movies.min_relatedness >= min_relatedness_sample_swap then true else false

    # merge aligned bams
    call Utils.merge_bams {
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
  }

  # somalier_relate_samples
  call Somalier.somalier_relate_samples {
    input:
      cohort_id = cohort_id,
      extracted_somalier_sites = flatten(somalier_extract.extracted_sites_per_sample),
      sample_ids = sample_id,
      coverages = mosdepth.mean_coverage,
      runtime_attributes = default_runtime_attributes
  }

  scatter (keep_drop in somalier_relate_samples.qc_related_keep_drop[1]) {
    Boolean pass_relatedness = if keep_drop == 'keep' then true else false
  }

  call Utils.summarize_qc {
    input:
      cohort_id = cohort_id,
      sample_ids = sample_id,
      n_movies = n_movies,
      qc_pass_sex = pass_sex,
      qc_pass_swap = pass_swap,
      qc_pass_relatedness = pass_relatedness,
      sex = somalier_relate_movies.inferred_sex,
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
    Array[IndexData] pbmm2_merged_bam = merged_bam
    File qc_summary = summarize_qc.quality_control_summary
    Array[Boolean] qc_pass_sex = pass_sex
    Array[Boolean] qc_pass_swap = pass_swap
    Array[Boolean] qc_pass_relatedness = pass_relatedness
  }
}
