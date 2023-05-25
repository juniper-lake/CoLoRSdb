version 1.0

import "../backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "../cohort_align_qc/cohort_align_qc.wdl" as CohortAlignQC
import "../cohort_align/cohort_align.wdl" as CohortAlign
import "../sample_call_variants/sample_call_variants.wdl" as SampleCallVariants
import "../cohort_combine_samples/cohort_combine_samples.wdl" as CohortCombineSamples

workflow cohort_main {
  input {
    Cohort cohort

    Array[ReferenceData] references
		Int primary_reference_idx = 0

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

	# align and QC samples to primary reference
	call CohortAlignQC.cohort_align_qc {
		input:
			cohort_id = cohort.cohort_id,
			samples = cohort.samples,
			max_sample_relatedness_qc = max_sample_relatedness_qc,
			min_movie_relatedness_qc = min_movie_relatedness_qc,
			reference = references[primary_reference_idx],
			default_runtime_attributes = default_runtime_attributes
	}

	# for each reference
	scatter (ref_idx in range(length(references))) {
		
		# sex and QC info come from the primary reference
		Array[AlignedSample] qc_samples = cohort_align_qc.aligned_samples
		
		# if not the primary reference, align the samples to the reference
		if (ref_idx != primary_reference_idx) {
			call CohortAlign.cohort_align {
				input:
					cohort_id = cohort.cohort_id,
					samples = cohort.samples,
					reference = references[ref_idx],
					default_runtime_attributes = default_runtime_attributes
			}
		}

		Array[AlignedSample] aligned_samples = select_first([cohort_align.aligned_samples, qc_samples])
		
		# for each sample
		scatter (sample_idx in range(length(qc_samples))) {

			# if the sample passed QC
			if (qc_samples[sample_idx].qc_pass) {

				String qc_pass_sample_id = qc_samples[sample_idx].sample_id
				String qc_pass_sex = qc_samples[sample_idx].sex
				IndexData qc_pass_aligned_bam = aligned_samples[sample_idx].aligned_bam
				
				# call variants for the sample
				call SampleCallVariants.sample_call_variants {
					input:
						sample_id = qc_pass_sample_id,
						aligned_bam = qc_pass_aligned_bam,
						sex = qc_pass_sex,
						reference = references[ref_idx],
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
				reference = references[ref_idx],
				default_runtime_attributes = default_runtime_attributes
		}
	}	

	output {
		File qc_summary_tsv = cohort_align_qc.qc_summary_tsv
		Array[File] coverage_summary_tsvs = select_all(cohort_align.coverage_summary_tsv)
		Array[IndexData] cohort_deepvariant_vcfs = cohort_combine_samples.cohort_deepvariant_vcf
		Array[IndexData] cohort_pbsv_vcfs = cohort_combine_samples.cohort_pbsv_vcf
		Array[IndexData] cohort_sniffles_vcfs = cohort_combine_samples.cohort_sniffles_vcf
		Array[IndexData] cohort_trgt_vcfs = cohort_combine_samples.cohort_trgt_vcf

		Array[File] pbsv_call_logs = flatten(cohort_combine_samples.pbsv_call_logs)
		Array[File?] hiphase_stats = cohort_combine_samples.hiphase_stats
		Array[File?] hiphase_blocks = cohort_combine_samples.hiphase_blocks
		Array[File?] hiphase_summaries = cohort_combine_samples.hiphase_summary
	}
}
