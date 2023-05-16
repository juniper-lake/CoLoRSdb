version 1.0

import "../backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "../sample_align_qc/sample_align_qc.wdl" as SampleAlignQC
import "../cohort_qc/cohort_qc.wdl" as CohortQC
# import "../sample_call_variants/sample_call_variants.wdl" as SampleCallVariants
# import "../cohort_combine_samples/cohort_combine_samples.wdl" as CohortCombineSamples

workflow cohort_main {
  input {
    Cohort cohort

    ReferenceData reference

    Float min_relatedness_sample_swap
		Float max_pairwise_relatedness

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


	scatter (sample in cohort.samples) {
		call SampleAlignQC.sample_align_qc {
			input:
				sample = sample,
				reference = reference,
				min_relatedness_sample_swap = min_relatedness_sample_swap,
				default_runtime_attributes = default_runtime_attributes
		}

		String sample_id = sample.sample_id
	}

	call CohortQC.cohort_qc {
		input:
			cohort_id = cohort.cohort_id,
			max_pairwise_relatedness = max_pairwise_relatedness,
			sample_ids = sample_id,
			extracted_somalier_sites = flatten(sample_align_qc.extracted_somalier_sites),
			n_movies = sample_align_qc.n_movies,
			qc_pass_sex = sample_align_qc.qc_pass_sex,
			qc_pass_swap = sample_align_qc.qc_pass_swap,
			sex = sample_align_qc.sex,
			min_movie_relatedness = sample_align_qc.min_movie_relatedness,
			coverage = sample_align_qc.coverage,
			read_count = sample_align_qc.read_count,
			unique_read_count = sample_align_qc.unique_read_count,
			read_quality_mean = sample_align_qc.read_quality_mean,
			read_quality_median = sample_align_qc.read_quality_median,
			read_quality_stdev = sample_align_qc.read_quality_stdev,
			read_length_mean = sample_align_qc.read_length_mean,
			read_length_median = sample_align_qc.read_length_median,
			read_length_stdev = sample_align_qc.read_length_stdev,
			default_runtime_attributes = default_runtime_attributes
	}

	# scatter (idx in range(length(cohort.samples))) {
	# 	if (cohort_qc.qc_pass) {
	# 		call SampleCallVariants.sample_call_variants {
	# 			input:

	# 		}
	# 	}
	# }

	# call CohortCombineSamples.cohort_combine_samples {
	# 	input:

	# }

	output {

	}
}

