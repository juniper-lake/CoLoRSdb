version 1.0

import "../structs.wdl"

task merge_bams {

  input {
    Array[File] bams
		String output_bam_name

		RuntimeAttributes runtime_attributes
  }

  Int threads = 8
	Int disk_size = ceil(size(bams[0], "GB") * length(bams) * 2 + 20)

  command {
		set -euo pipefail

		if [[ "~{length(bams)}" -eq 1 ]]; then
			mv ~{bams[0]} ~{output_bam_name}
		else
			samtools merge \
				-@ ~{threads - 1} \
				-o ~{output_bam_name} \
				~{sep=' ' bams}
		fi

		samtools index ~{output_bam_name}
  }

  output {
    File merged_bam = output_bam_name
    File merged_bam_index = "~{output_bam_name}.bai"
  }

  runtime {
    cpu: threads
    memory: "1 GB"
		disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/samtools:1.14"
  }
}


task zip_index_vcf {

  input {
    File vcf

		RuntimeAttributes runtime_attributes
  }

	Int threads = 4
	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)
  String vcf_basename = basename(vcf)

  command <<<
    set -euo pipefail

    bgzip \
			--threads ~{threads} \
			~{vcf} -c \
			> ~{vcf_basename}.gz

    tabix \
			--preset vcf \
			~{vcf_basename}.gz
  >>>

  output {
    File zipped_vcf = "~{vcf_basename}.gz"
    File zipped_vcf_index = "~{vcf_basename}.gz.tbi"
  }

  runtime {
    cpu: threads
    memory: "4 GB"
		disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/htslib:1.14"
  }
}


task index_vcf {

  input {
    File gzipped_vcf

		RuntimeAttributes runtime_attributes
  }

	Int threads = 1
	Int disk_size = ceil(size(gzipped_vcf, "GB") * 2 + 20)
  String vcf_basename = basename(gzipped_vcf)

  command <<<
    set -euo pipefail

    tabix \
			--preset vcf \
			~{vcf_basename}.gz
  >>>

  output {
    File zipped_vcf = "~{vcf_basename}.gz"
    File zipped_vcf_index = "~{vcf_basename}.gz.tbi"
  }

  runtime {
    cpu: threads
    memory: "2 GB"
		disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/htslib:1.14"
  }
}
