version 1.0

import "backend_configuration.wdl" as BackendConfiguration
import "align_qc.wdl" as AlignQC
import "call_variants_by_sample.wdl" as CallVariantsBySample
import "merge_samples.wdl" as MergeSamples
import "../tasks/utils.wdl" as Utils

workflow colors_main {
  input {
    String cohort_id
    File sample_sheet
    Boolean align_qc_only = false

    File reference_bundle
    Boolean anonymize_output = false

    # backend configuration
    String backend
    String? zones
    String? aws_spot_queue_arn
    String? aws_on_demand_queue_arn
    Boolean preemptible
  }

  # these inputs are not in the input section because should not be changed for data contributing to CoLoRSdb
  String container_registry = "quay.io/colorsdb"

  call BackendConfiguration.backend_configuration {
    input:
      backend = backend,
      zones = zones,
      aws_spot_queue_arn = aws_spot_queue_arn,
      aws_on_demand_queue_arn = aws_on_demand_queue_arn,
      container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

  # unzip reference bundle
  call Utils.unzip_reference_bundle {
    input:
      reference_bundle = reference_bundle,
      runtime_attributes = default_runtime_attributes
  }

  call Utils.read_sample_sheet {
    input:
      sample_sheet = sample_sheet,
      runtime_attributes = default_runtime_attributes
  }

  scatter (sample_idx in range(length(read_sample_sheet.sample_ids))) {
    # override sample qc as a determinant of if variants are called
    if (read_sample_sheet.qc_pass[sample_idx]!="null") {
      Boolean override_qc_pass = if read_sample_sheet.qc_pass[sample_idx] == "true" then true else false
    }
    # override sex
    if (read_sample_sheet.sexes[sample_idx] != "null") {
      String override_sex = read_sample_sheet.sexes[sample_idx]
    }
  }

  scatter (reference in [unzip_reference_bundle.grch38, unzip_reference_bundle.chm13]) {
    # align and qc
    call AlignQC.align_qc {
      input:
        cohort_id = cohort_id,
        sample_ids = read_sample_sheet.sample_ids,
        sample_movies = read_sample_sheet.movies,
        reference = reference,
        default_runtime_attributes = default_runtime_attributes
    }

    if (align_qc_only) {
      Array[File] output_aligned_bams = align_qc.aligned_bams
      Array[File] output_aligned_bam_indexes = align_qc.aligned_bam_indexes
    }

    scatter (idx in range(length(read_sample_sheet.sample_ids))) {
      if (select_first([override_qc_pass[idx], align_qc.qc_pass[idx]])) {
        String qc_pass_sample_id = read_sample_sheet.sample_ids[idx]
        String qc_pass_sex = select_first([override_sex[idx], align_qc.sexes[idx]])
        File qc_pass_bam = align_qc.aligned_bams[idx]
        File qc_pass_bam_index = align_qc.aligned_bam_indexes[idx]
      }
    }

    Int n_samples = length(select_all(qc_pass_sample_id))

    if (!align_qc_only) {
      # continue only if more than one sample had variants called
      if (n_samples > 1) {
        scatter (idx in range(n_samples)) {
          call CallVariantsBySample.call_variants_by_sample {
            input:
              sample_id = select_all(qc_pass_sample_id)[idx],
              sex = select_all(qc_pass_sex)[idx],
              aligned_bam = select_all(qc_pass_bam)[idx],
              aligned_bam_index = select_all(qc_pass_bam_index)[idx],
              reference = reference,
              default_runtime_attributes = default_runtime_attributes
          }
        }

        # merge samples
        call MergeSamples.merge_samples {
          input:
            cohort_id = cohort_id,
            anonymize_output = anonymize_output,
            sample_ids = select_all(qc_pass_sample_id),
            sexes = select_all(qc_pass_sex),
            pbsv_vcfs = select_all(call_variants_by_sample.pbsv_vcf),
            pbsv_vcf_indexes = select_all(call_variants_by_sample.pbsv_vcf_index),
            unzipped_pbsv_vcfs = select_all(call_variants_by_sample.unzipped_pbsv_vcf),
            deepvariant_gvcfs = select_all(call_variants_by_sample.deepvariant_gvcf),
            deepvariant_gvcf_indexes = select_all(call_variants_by_sample.deepvariant_gvcf_index),
            sniffles_snfs = select_all(call_variants_by_sample.sniffles_snf),
            trgt_vcfs = select_all(call_variants_by_sample.trgt_vcf),
            reference = reference,
            default_runtime_attributes = default_runtime_attributes
        }
      }
    }
  }

  output {
    # messages
    String message = if (select_first(n_samples) < 2) then "At least two samples passing QC are required to call variants." else "Variants were called and merged for ~{select_first(n_samples)} samples."

    # align and qc
    Array[File] somalier_pairs = align_qc.somalier_pairs
    Array[File] somalier_samples = align_qc.somalier_samples
    Array[File] sample_qc_summary_tsv = align_qc.sample_qc_summary_tsv
    Array[File] bam_qc_summary_tsv = align_qc.bam_qc_summary_tsv
    Array[Array[File]?]+ aligned_bams = output_aligned_bams
    Array[Array[File]?]+ aligned_bam_indexes = output_aligned_bam_indexes

    # postprocessed VCFs
    Array[File?]+ deepvariant_glnexus_postprocessed_vcf = merge_samples.deepvariant_glnexus_postprocessed_vcf
    Array[File?]+ deepvariant_glnexus_postprocessed_vcf_index = merge_samples.deepvariant_glnexus_postprocessed_vcf_index
    Array[File?]+ pbsv_jasminesv_postprocessed_vcf = merge_samples.pbsv_jasminesv_postprocessed_vcf
    Array[File?]+ pbsv_jasminesv_postprocessed_vcf_index = merge_samples.pbsv_jasminesv_postprocessed_vcf_index
    Array[File?]+ sniffles_postprocessed_vcf = merge_samples.sniffles_postprocessed_vcf
    Array[File?]+ sniffles_postprocessed_vcf_index = merge_samples.sniffles_postprocessed_vcf_index
    Array[Array[File]?]+ trgt_postprocessed_vcfs = merge_samples.trgt_postprocessed_vcfs
    Array[Array[File]?]+ trgt_postprocessed_vcf_indexes = merge_samples.trgt_postprocessed_vcf_indexes

    # original VCFs, so if anonymize_output=false we can access VCFs without ploidy changes
    Array[File?]+ deepvariant_glnexus_vcf = merge_samples.deepvariant_glnexus_vcf
    Array[File?]+ deepvariant_glnexus_vcf_index = merge_samples.deepvariant_glnexus_vcf_index
    Array[File?]+ sniffles_vcf = merge_samples.sniffles_vcf
    Array[File?]+ sniffles_vcf_index = merge_samples.sniffles_vcf_index
    Array[File?]+ pbsv_jasminesv_vcf = merge_samples.pbsv_jasminesv_vcf
    Array[File?]+ pbsv_jasminesv_vcf_index = merge_samples.pbsv_jasminesv_vcf_index
    Array[Array[File]?]+ trgt_vcf = merge_samples.trgt_vcf
    Array[Array[File]?]+ trgt_vcf_index = merge_samples.trgt_vcf_index

    # vcf stats
    Array[File?]+ pbsv_jasminesv_vcf_stats = merge_samples.pbsv_jasminesv_vcf_stats
    Array[File?]+ sniffles_vcf_stats = merge_samples.sniffles_vcf_stats

    # ancestry when peddy is run
    Array[File?]+ peddy_het_check = merge_samples.peddy_het_check
    Array[File?]+ peddy_sex_check = merge_samples.peddy_sex_check
    Array[File?]+ peddy_ped_check = merge_samples.peddy_ped_check
    Array[File?]+ peddy_background_pca = merge_samples.peddy_background_pca
    Array[File?]+ peddy_html = merge_samples.peddy_html
    Array[File?]+ peddy_ped = merge_samples.peddy_ped
    Array[File?]+ peddy_vs_html = merge_samples.peddy_vs_html
  }
}
