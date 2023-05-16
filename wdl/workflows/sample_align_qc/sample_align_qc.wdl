version 1.0

import "../../tasks/pbmm2.wdl" as Pbmm2
import "../../tasks/somalier.wdl" as Somalier
import "../../tasks/utils.wdl" as Utils
import "../../tasks/mosdepth.wdl" as Mosdepth


workflow sample_align_qc {

  input {
    Sample sample

    ReferenceData reference
    Float min_relatedness_sample_swap

    RuntimeAttributes default_runtime_attributes
  }

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

  # check for sample swaps among movies and confirm sex of sample
  call Somalier.somalier_relate_movies {
    input:
      sample_id = sample.sample_id,
      extracted_somalier_sites = somalier_extract.extracted_sites_per_movie,
      runtime_attributes = default_runtime_attributes
  }

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


  output {
    IndexData merged_aligned_bam = merged_bam
    String sex = somalier_relate_movies.inferred_sex
    Float min_movie_relatedness = somalier_relate_movies.min_relatedness
    Boolean qc_pass_sex = if somalier_relate_movies.inferred_sex != "UNKNOWN" then true else false
    Boolean qc_pass_swap = if somalier_relate_movies.min_relatedness >= min_relatedness_sample_swap then true else false
    Int n_movies = length(sample.movies)
    Float coverage = mosdepth.mean_coverage
    Array[File] extracted_somalier_sites = somalier_extract.extracted_sites_per_sample
    Int read_count = combine_smrtcell_stats.read_count
    Int unique_read_count = combine_smrtcell_stats.unique_read_count
    Float read_quality_mean = combine_smrtcell_stats.read_quality_mean
    Float read_quality_median = combine_smrtcell_stats.read_quality_median
    Float read_quality_stdev = combine_smrtcell_stats.read_quality_stdev
    Float read_length_mean = combine_smrtcell_stats.read_length_mean
    Float read_length_median = combine_smrtcell_stats.read_length_median
    Float read_length_stdev = combine_smrtcell_stats.read_length_stdev
  }
}
