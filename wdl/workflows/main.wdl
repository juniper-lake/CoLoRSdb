version 1.0

import "./backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "../tasks/utils.wdl" as Utils
import "./cohort_align_qc/cohort_align_qc.wdl" as CohortAlignQC
import "./sample_call_variants/sample_call_variants.wdl" as SampleCallVariants
import "./cohort_combine_samples/cohort_combine_samples.wdl" as CohortCombineSamples

workflow colors_cohort {
  input {
    String cohort_id
    File sample_sheet

    ReferenceData reference
    Boolean anonymize_output = true

    # backend configuration
    String backend
    String? zones
    String? aws_spot_queue_arn
    String? aws_on_demand_queue_arn
    Boolean preemptible
  }

  # these inputs are not in the input section because should not be changed for data contributing to CoLoRSdb
  String container_registry = "quay.io/colorsdb"
  Float max_sample_relatedness_qc = 0.125
  Float min_movie_relatedness_qc = 0.7
  String deepvariant_version = "1.5.0"

  call BackendConfiguration.backend_configuration {
    input:
      backend = backend,
      zones = zones,
      aws_spot_queue_arn = aws_spot_queue_arn,
      aws_on_demand_queue_arn = aws_on_demand_queue_arn,
      container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

  # read sample sheet
  call Utils.read_sample_sheet {
    input:
      sample_sheet = sample_sheet,
      runtime_attributes = default_runtime_attributes
  }

  scatter (sample in read_sample_sheet.samples) {
    Sample samples = object {
      sample_id: sample.left,
      movies: sample.right
    }
  }

  # require at least 2 samples to run alignment and QC
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
  
  Array[AlignedSample] qc_pass_samples = select_first([cohort_align_qc.qc_pass_samples,[]])

  # continue only if more than one sample passes QC
  if (length(qc_pass_samples) > 1) {
    # for each sample
    scatter (qc_pass_sample in qc_pass_samples) {
      String sample_id = qc_pass_sample.sample_id
      String sex = qc_pass_sample.sex
      # call variants for the sample
      call SampleCallVariants.sample_call_variants {
        input:
          sample_id = sample_id,
          sex = sex,
          aligned_bam = qc_pass_sample.aligned_bam,
          reference = reference,
          deepvariant_version = deepvariant_version,
          default_runtime_attributes = default_runtime_attributes
      }
    }

    # combine variant calls from samples that passed QC
    call CohortCombineSamples.cohort_combine_samples {
      input:
        cohort_id = cohort_id,
        anonymize_output = anonymize_output,
        sample_ids = sample_id,
        sexes = sex,
        pbsv_vcfs = sample_call_variants.pbsv_vcf,
        unzipped_pbsv_vcfs = sample_call_variants.unzipped_pbsv_vcf,
        deepvariant_gvcfs = sample_call_variants.deepvariant_gvcf,
        sniffles_snfs = sample_call_variants.sniffles_snf,
        trgt_vcfs = select_all(sample_call_variants.trgt_vcf),
        hificnv_vcfs = select_all(sample_call_variants.hificnv_vcf),
        reference = reference,
        default_runtime_attributes = default_runtime_attributes
    }
  }

  output {
    # messages
    String sample_size_message = if (length(samples) < 2) then "At least two samples are required to run the workflow, but only ~{length(samples)} samples were included in cohort." else "~{length(samples)} samples are included in alignment and QC steps."
    String qc_message = if (length(qc_pass_samples) < 2) then "At least two samples are required to pass QC, but only ~{length(qc_pass_samples)} samples passed QC." else "~{length(qc_pass_samples)} samples out of ~{length(samples)} samples passed QC and will be used for variant calling."

    # qc
    File? qc_summary_tsv = cohort_align_qc.qc_summary_tsv
    File? somalier_pairs_tsv = cohort_align_qc.somalier_pairs
    File? somalier_samples_tsv = cohort_align_qc.somalier_samples

    # postprocessed vcfs
    IndexData? deepvariant_glnexus_postprocessed_vcf = cohort_combine_samples.deepvariant_glnexus_postprocessed_vcf
    IndexData? pbsv_jasminesv_postprocessed_vcf = cohort_combine_samples.pbsv_jasminesv_postprocessed_vcf
    IndexData? sniffles_postprocessed_vcf = cohort_combine_samples.sniffles_postprocessed_vcf
    Array[IndexData]? trgt_postprocessed_vcf = cohort_combine_samples.trgt_postprocessed_vcf
    IndexData? hificnv_postprocessed_vcf = cohort_combine_samples.hificnv_postprocessed_vcf

    # original vcfs
    IndexData? deepvariant_glnexus_vcf = cohort_combine_samples.deepvariant_glnexus_vcf
    IndexData? pbsv_strictmerge_vcf = cohort_combine_samples.pbsv_strictmerge_vcf
    IndexData? pbsv_jasminesv_vcf = cohort_combine_samples.pbsv_jasminesv_vcf
    IndexData? sniffles_vcf = cohort_combine_samples.sniffles_vcf
    Array[IndexData]? trgt_vcf = cohort_combine_samples.trgt_vcf
    IndexData? hificnv_vcf = cohort_combine_samples.hificnv_vcf

    # vcf stats
    File? pbsv_vcf_stats = cohort_combine_samples.pbsv_jasminesv_vcf_stats
    File? sniffles_vcf_stats = cohort_combine_samples.sniffles_vcf_stats
    
    # ancestry
    File? peddy_het_check = cohort_combine_samples.peddy_het_check
    File? peddy_sex_check = cohort_combine_samples.peddy_sex_check
    File? peddy_ped_check = cohort_combine_samples.peddy_ped_check
    File? peddy_background_pca = cohort_combine_samples.peddy_background_pca
    File? peddy_html = cohort_combine_samples.peddy_html
    File? peddy_ped = cohort_combine_samples.peddy_ped
    File? peddy_vs_html = cohort_combine_samples.peddy_vs_html
  }
}
