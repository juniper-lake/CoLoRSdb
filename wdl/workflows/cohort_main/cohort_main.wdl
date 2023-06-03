version 1.0

import "../backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "../cohort_align_qc/cohort_align_qc.wdl" as CohortAlignQC
import "../sample_call_variants/sample_call_variants.wdl" as SampleCallVariants
import "../cohort_combine_samples/cohort_combine_samples.wdl" as CohortCombineSamples

workflow cohort_main {
  input {
    Cohort cohort

    ReferenceData reference

    Float max_sample_relatedness_qc = 0.125
    Float min_movie_relatedness_qc = 0.875

		String deepvariant_version

		# Backend configuration
		String backend
		String? zones
		String? aws_spot_queue_arn
		String? aws_on_demand_queue_arn

		Boolean preemptible

  }

	call BackendConfiguration.backend_configuration {
		input:
			backend = backend,
			zones = zones,
			aws_spot_queue_arn = aws_spot_queue_arn,
			aws_on_demand_queue_arn = aws_on_demand_queue_arn
	}

	RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

	call CohortAlignQC.cohort_align_qc {
		input:
			cohort_id = cohort.cohort_id,
			samples = cohort.samples,
			max_sample_relatedness_qc = max_sample_relatedness_qc,
			min_movie_relatedness_qc = min_movie_relatedness_qc,
			reference = reference,
			default_runtime_attributes = default_runtime_attributes
	}
	
	# for each sample
	scatter (sample_idx in range(length(cohort.samples))) {

		# if the sample passed QC
		if (cohort_align_qc.qc_pass[sample_idx]) {
			
			String qc_pass_sample_id = cohort.samples[sample_idx].sample_id
			IndexData qc_pass_aligned_bam = cohort_align_qc.aligned_bams[sample_idx]
			String qc_pass_sex = cohort_align_qc.inferred_sexes[sample_idx]

			# call variants for the sample
			call SampleCallVariants.sample_call_variants {
				input:
					sample_id = qc_pass_sample_id,
					sex = qc_pass_sex,
					aligned_bam = qc_pass_aligned_bam,
					reference = reference,
					deepvariant_version = deepvariant_version,
					default_runtime_attributes = default_runtime_attributes
			}
		}
	}

	# combine variant calls from samples that passed QC
	call CohortCombineSamples.cohort_combine_samples {
		input:
			cohort_id = cohort.cohort_id,
			anonymize_output = cohort.anonymize_output,
			sample_ids = select_all(qc_pass_sample_id),
			sexes = select_all(qc_pass_sex),
			aligned_bams = select_all(qc_pass_aligned_bam),
			svsigs = select_all(sample_call_variants.pbsv_svsigs),
			gvcfs = select_all(sample_call_variants.deepvariant_gvcf),
			snfs = select_all(sample_call_variants.sniffles_snf),
			trgt_vcfs = select_all(sample_call_variants.trgt_vcf),
			reference = reference,
			default_runtime_attributes = default_runtime_attributes
	}

	output {
		File qc_summary_tsv = cohort_align_qc.qc_summary_tsv
		File pairwise_relatedness_tsv = cohort_align_qc.pairwise_relatedness

		IndexData cohort_deepvariant_vcf = cohort_combine_samples.cohort_deepvariant_vcf
		IndexData cohort_pbsv_vcf = cohort_combine_samples.cohort_pbsv_vcf
		IndexData cohort_sniffles_vcf = cohort_combine_samples.cohort_sniffles_vcf
		IndexData? cohort_trgt_vcf = cohort_combine_samples.cohort_trgt_vcf
		IndexData? cohort_hificnv_vcf = cohort_combine_samples.cohort_hificnv_vcf

		File? peddy_het_check = cohort_combine_samples.peddy_het_check
    File? peddy_sex_check = cohort_combine_samples.peddy_sex_check
    File? peddy_ped_check = cohort_combine_samples.peddy_ped_check
    File? peddy_background_pca = cohort_combine_samples.peddy_background_pca
    File? peddy_html = cohort_combine_samples.peddy_html
    File? peddy_ped = cohort_combine_samples.peddy_ped
    File? peddy_vs_html = cohort_combine_samples.peddy_vs_html

		Array[File] pbsv_call_logs = cohort_combine_samples.pbsv_call_logs
		File? hiphase_stats = cohort_combine_samples.hiphase_stats
		File? hiphase_blocks = cohort_combine_samples.hiphase_blocks
		File? hiphase_summary = cohort_combine_samples.hiphase_summary
	}
}
