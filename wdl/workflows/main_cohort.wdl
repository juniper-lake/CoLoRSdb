version 1.0

import "./backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "../tasks/utils.wdl" as Utils
import "./cohort_align_qc/cohort_align_qc.wdl" as CohortAlignQC
import "./sample_call_variants/sample_call_variants.wdl" as SampleCallVariants
import "./cohort_combine_samples/cohort_combine_samples.wdl" as CohortCombineSamples

workflow cohort_main {
  input {
    String cohort_id
    Array[Sample]? cohort_samples
    File? cohort_sample_sheet

    ReferenceData reference

    Boolean anonymize_output = true
    Int max_samples_pbsv_call = 150
    Float max_sample_relatedness_qc = 0.125
    Float min_movie_relatedness_qc = 0.875

    String deepvariant_version

    # Backend configuration
    String backend
    String? zones
    String? aws_spot_queue_arn
    String? aws_on_demand_queue_arn
    Boolean preemptible
    String container_registry
  }

  call BackendConfiguration.backend_configuration {
    input:
      backend = backend,
      zones = zones,
      aws_spot_queue_arn = aws_spot_queue_arn,
      aws_on_demand_queue_arn = aws_on_demand_queue_arn,
      container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

  if (!defined(cohort_samples)) {
    call Utils.read_sample_sheet {
      input:
        sample_sheet = select_first([cohort_sample_sheet]),
        runtime_attributes = default_runtime_attributes
    }

    scatter (sample_idx in range(length(read_sample_sheet.sample_ids))) {
      Sample sample_sheet_samples = object {
        "sample_id": read_sample_sheet.sample_ids[sample_idx],
        "movies": read_sample_sheet.movies[sample_idx]
      }
    }
  }
  
  Array[Sample] samples = select_first([cohort_samples, sample_sheet_samples])

  if (length(samples) > 1) {
    call CohortAlignQC.cohort_align_qc {
    input:
      cohort_id = cohort_id,
      samples = samples,
      max_sample_relatedness_qc = max_sample_relatedness_qc,
      min_movie_relatedness_qc = min_movie_relatedness_qc,
      reference = reference,
      default_runtime_attributes = default_runtime_attributes
    }
  }

  # continue only if more than one sample passes QC
  if (length(select_first([cohort_align_qc.qc_pass_sample_ids, []])) > 1) {
    # for each sample
    scatter (sample_idx in range(length(samples))) {
      # if the sample passed QC
      if (cohort_align_qc.qc_pass[sample_idx]) {
        
        String sample_id = samples[sample_idx].sample_id
        IndexData aligned_bam = cohort_align_qc.aligned_bams[sample_idx]
        String sex = cohort_align_qc.inferred_sexes[sample_idx]

        # call variants for the sample
        call SampleCallVariants.sample_call_variants {
          input:
            sample_id = sample_id,
            sex = sex,
            aligned_bam = aligned_bam,
            reference = reference,
            deepvariant_version = deepvariant_version,
            default_runtime_attributes = default_runtime_attributes
        }
      }
    }

    # combine variant calls from samples that passed QC
    call CohortCombineSamples.cohort_combine_samples {
      input:
        cohort_id = cohort_id,
        anonymize_output = anonymize_output,
        max_samples_pbsv_call = max_samples_pbsv_call,
        sample_ids = select_all(sample_id),
        sexes = select_all(sex),
        aligned_bams = select_all(aligned_bam),
        pbsv_svsigs = select_all(sample_call_variants.pbsv_svsigs),
        deepvariant_gvcfs = select_all(sample_call_variants.deepvariant_gvcf),
        sniffles_snfs = select_all(sample_call_variants.sniffles_snf),
        trgt_vcfs = select_all(sample_call_variants.trgt_vcf),
        hificnv_vcfs = select_all(sample_call_variants.hificnv_vcf),
        reference = reference,
        default_runtime_attributes = default_runtime_attributes
    }
  }

  output {
    File? qc_summary_tsv = cohort_align_qc.qc_summary_tsv
    File? pairwise_relatedness_tsv = cohort_align_qc.pairwise_relatedness

    IndexData? cohort_deepvariant_vcf = cohort_combine_samples.cohort_deepvariant_vcf
    Array[IndexData]? cohort_pbsv_vcf = cohort_combine_samples.cohort_pbsv_vcfs
    IndexData? cohort_sniffles_vcf = cohort_combine_samples.cohort_sniffles_vcf
    IndexData? cohort_trgt_vcf = cohort_combine_samples.cohort_trgt_vcf
    IndexData? cohort_hificnv_vcf = cohort_combine_samples.cohort_hificnv_vcf

    File? cohort_deepvariant_vcf_stats = cohort_combine_samples.cohort_deepvariant_vcf_stats
    Array[File]? cohort_pbsv_vcf_stats = cohort_combine_samples.cohort_pbsv_vcf_stats
    File? cohort_sniffles_vcf_stats = cohort_combine_samples.cohort_sniffles_vcf_stats

    File? peddy_het_check = cohort_combine_samples.peddy_het_check
    File? peddy_sex_check = cohort_combine_samples.peddy_sex_check
    File? peddy_ped_check = cohort_combine_samples.peddy_ped_check
    File? peddy_background_pca = cohort_combine_samples.peddy_background_pca
    File? peddy_html = cohort_combine_samples.peddy_html
    File? peddy_ped = cohort_combine_samples.peddy_ped
    File? peddy_vs_html = cohort_combine_samples.peddy_vs_html
  }
}
