version 1.0

import "../../tasks/pbmm2.wdl" as Pbmm2
import "../../tasks/somalier.wdl" as Somalier
import "../../tasks/utils.wdl" as Utils


workflow align_and_qc {

  input {
    String cohort_id
    Array[Sample] samples

    ReferenceData reference
    Float min_relatedness_sample_swap

    RuntimeAttributes default_runtime_attributes
  }

  scatter (sample in samples) {

    # align each movie with pbmm2
    scatter (idx in range(length(sample.movies))) {     
      call Pbmm2.pbmm2_align {
          input: 
            movie = sample.movies[idx],
            sample_id = sample.sample_id,
            reference_name = reference.name,
            reference_fasta = reference.fasta.data,
            reference_index = reference.fasta.index,
            runtime_attributes = default_runtime_attributes
      }
      
      IndexData aligned_bam = {
        "data": pbmm2_align.bam,
        "index": pbmm2_align.bam_index
      }

      # get movie name from movie path
      String movie_name = sub(basename(sample.movies[idx]), "\\..*", "")

      # extract somalier sites from aligned bam
      call Somalier.somalier_extract {
        input:
          sample_id = sample.sample_id,
          movie_name = movie_name,
          bam = aligned_bam.data,
          bam_index = aligned_bam.index,
          reference_fasta = reference.fasta.data,
          reference_index = reference.fasta.index,
          somalier_sites_vcf = reference.somalier_sites_vcf,
          runtime_attributes = default_runtime_attributes
      }
    }

    # check for sample swaps and confirm gender of sample
    call Somalier.somalier_relate_movies {
      input:
        sample_id = sample.sample_id,
        extracted_somalier_sites = somalier_extract.extracted_sites,
        runtime_attributes = default_runtime_attributes
    }

    Boolean pass_sex = if somalier_relate_movies.inferred_sex != "UNKNOWN" then true else false
    Boolean pass_swap = if somalier_relate_movies.min_relatedness >= min_relatedness_sample_swap then true else false

    if (somalier_relate_movies.inferred_sex != "UNKNOWN" && somalier_relate_movies.min_relatedness >= min_relatedness_sample_swap) {
      # merge aligned bams for trgt and pbsv
      call Utils.merge_bams {
        input:
          bams = pbmm2_align.bam,
          output_bam_name = "~{sample.sample_id}.~{reference.name}.bam",
          runtime_attributes = default_runtime_attributes
      }

      IndexData merged_bam = {
        "data": merge_bams.merged_bam, 
        "index": merge_bams.merged_bam_index
      }
    }

    # run somalier extract on merged bam
    # run somalier relate with script to check for related samples
    # run script to identify which samples need to be removed
    # run mosdepth on merged bam -> good to know coverage? what output? just one number?


  }

  output {
    # QC pass
    Array[IndexData] pbmm2_merged_bam = select_all(merged_bam)
    # File qc_log = quality_control_log.summary

    # QC fail
    Array[Boolean] qc_pass_sex = pass_sex
    Array[Boolean] qc_pass_swap = pass_swap
  }
}
