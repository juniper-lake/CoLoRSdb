version 1.0

import "../structs.wdl"

task hificnv {
	input {
		String sample_id
		String sex

		File bam
		File bam_index

		File small_variant_vcf
		File small_variant_vcf_index

    String reference_name
		File reference
		File reference_index

		File exclude_bed
		File exclude_bed_index

		File expected_bed_male
		File expected_bed_female

		RuntimeAttributes runtime_attributes
	}

  String output_prefix = "~{sample_id}.~{reference_name}.hificnv"

	Boolean sex_defined = defined(sex)
	File expected_bed = if sex == "male" then expected_bed_male else expected_bed_female

	Int threads = 8
	Int mem_gb = threads * 2
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB"))+ 20)

	command <<<
		set -euo pipefail

		hificnv --version

		hificnv \
			--threads ~{threads} \
			--bam ~{bam} \
			--ref ~{reference} \
			--maf ~{small_variant_vcf} \
			--exclude ~{exclude_bed} \
			--expected-cn ~{expected_bed} \
			--output-prefix ~{output_prefix}
	>>>

	output {
		File cnv_vcf = "~{output_prefix}.~{sample_id}.vcf.gz"
	}

	runtime {
		cpu: threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
		docker: "~{runtime_attributes.container_registry}/hificnv:0.1.6"
	}
}
