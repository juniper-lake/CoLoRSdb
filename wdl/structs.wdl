version 1.0

struct IndexData {
	File data
	File index
}

struct Sample {
	String sample_id
	Array[File] movies
}

struct AlignedSample {
	String sample_id
	IndexData aligned_bam
	String sex
	Boolean qc_pass
}

struct Cohort {
	String cohort_id
	Array[Sample] samples

	Boolean anonymize_output
}

struct ReferenceData {
	String name
	IndexData fasta

	Array[String] chromosomes

	File tandem_repeat_bed
	File trgt_tandem_repeat_bed

	File somalier_sites_vcf
	File? peddy_sites
	File? peddy_bin
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
