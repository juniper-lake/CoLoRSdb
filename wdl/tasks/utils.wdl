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
    memory: "1GB"
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