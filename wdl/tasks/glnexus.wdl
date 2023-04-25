version 1.0

import "../structs.wdl"

task glnexus {
	input {
		String cohort_id
		Array[File] gvcfs
		Array[File] gvcf_indexes

		String reference_name

		File? regions_bed

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gb = 30
	Int disk_size = ceil((size(gvcfs[0], "GB") * length(gvcfs)) * 2 + 100)

  String output_prefix =  "~{cohort_id}.~{reference_name}.deepvariant.glnexus"

	command <<<
		set -euo pipefail

		glnexus_cli \
			--threads ~{threads} \
			--mem-gbytes ~{mem_gb} \
			--dir ~{cohort_id}.~{reference_name}.GLnexus.DB \
			--config DeepVariant_unfiltered \
			~{"--bed " + regions_bed} \
			~{sep=' ' gvcfs} \
		> ~{output_prefix}.bcf

		bcftools view \
			--threads ~{threads} \
			--output-type z \
			--output-file ~{output_prefix}.vcf.gz \
			~{output_prefix}.bcf

		tabix ~{output_prefix}.vcf.gz
	>>>

	output {
		File vcf = "~{output_prefix}.vcf.gz"
		File vcf_index = "~{output_prefix}.vcf.gz.tbi"
	}

	runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
		disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
		docker: "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
	}
}