version 1.0

struct IndexData {
	File data
	File index
}

struct Sample {
	String sample_id
	Array[File] movies
}

struct Cohort {
	String cohort_id
	Array[Sample] samples

	Boolean aggregate_output
}

struct ReferenceData {
	String name
	IndexData fasta

	Array[String] chromosomes
	File chromosome_lengths

	File tandem_repeat_bed
	File trgt_tandem_repeat_bed

	File somalier_sites_vcf
}

struct RuntimeAttributes {
	# The number of times to retry a task that fails due to preemption
	Int preemptible_tries
	# The number of times to retry a task that fails due a to nonzero return code
	Int max_retries

	String zones
	String queue_arn
	String container_registry
}
