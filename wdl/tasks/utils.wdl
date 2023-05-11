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
    cpu: 1
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


task concat_vcfs {
  
  input {
    Array[File] vcfs
    String output_vcf_name

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs[0], "GB") * length(vcfs) * 2 + 20)

  command {
    set -euo pipefail

    if [[ "~{length(vcfs)}" -eq 1 ]]; then
      mv ~{vcfs[0]} ~{output_vcf_name}
    else
      # will throw error if vcfs are not in correct order
      # could --allow-overlaps but would require sort after and that's compute intensive
      bcftools concat \
        -o ~{output_vcf_name} \
        ~{sep=' ' vcfs}
    fi
  }

  output {
    File concatenated_vcf = output_vcf_name
  }

  runtime {
    cpu: 1
    memory: "2 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/bcftools:1.14"
  }
}


task combine_smrtcell_stats {
  input {
    String sample_id
    Array[File] read_length_and_quality_tsvs

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(read_length_and_quality_tsvs, "GB") + 20)

  command <<<
    set -euo pipefail

    cat ~{sep=' ' read_length_and_quality_tsvs} \
      | datamash count 1 countunique 1 median 2 mean 2 pstdev 2 median 3 mean 3 pstdev 3 \
      | awk '{OFS="\n"; $1=$1}1'
      > ~{sample_id}.smrtcell_stats.txt
    
    sed '1q;d' ~{sample_id}.smrtcell_stats.txt > read_count.txt
    sed '2q;d' ~{sample_id}.smrtcell_stats.txt > read_count_unique.txt
    sed '3q;d' ~{sample_id}.smrtcell_stats.txt > read_length_median.txt
    sed '4q;d' ~{sample_id}.smrtcell_stats.txt > read_length_mean.txt
    sed '5q;d' ~{sample_id}.smrtcell_stats.txt > read_length_stdev.txt
    sed '6q;d' ~{sample_id}.smrtcell_stats.txt > read_quality_median.txt
    sed '7q;d' ~{sample_id}.smrtcell_stats.txt > read_quality_mean.txt
    sed '8q;d' ~{sample_id}.smrtcell_stats.txt > read_quality_stdev.txt

  >>>

  output {
    Int read_count = read_int("read_count.txt")
    Int unique_read_count = read_int("read_count_unique.txt")
    Float read_length_median = read_float("read_length_median.txt")
    Float read_length_mean = read_float("read_length_mean.txt")
    Float read_length_stdev = read_float("read_length_stdev.txt")
    Float read_quality_median = read_float("read_quality_median.txt")
    Float read_quality_mean = read_float("read_quality_mean.txt")
    Float read_quality_stdev = read_float("read_quality_stdev.txt")
  }

  runtime {
    cpu: 1
    memory: "2 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "datamash:1.1.0"
  }  
}


task summarize_qc {
  input {
    String cohort_id
    Array[String] sample_ids
    Array[Int] n_movies
    Array[Boolean] qc_pass_sex
    Array[Boolean] qc_pass_swap
    Array[Boolean] qc_pass_relatedness
    Array[String] sex
    Array[Float] coverage_mean
    Array[Float] read_count
    Array[Float] unique_read_count
    Array[Float] read_quality_mean
    Array[Float] read_quality_median
    Array[Float] read_quality_stdev
    Array[Float] read_length_mean
    Array[Float] read_length_median
    Array[Float] read_length_stdev

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 20

  command <<<
    set -euo pipefail
    echo -e "sample_id\tn_movies\tqc_pass_sex\tqc_pass_swap\tqc_pass_relatedness\tsex\tcoverage_mean\tread_count\tunique_read_count\tread_quality_mean\tread_quality_median\tread_quality_stdev\tread_length_mean\tread_length_median\tread_length_stdev" \
      > ~{cohort_id}.quality_control_summary.tsv
    
    paste <(echo -e ~{sep="\n" sample_ids}) \
      <(echo -e ~{sep="\n" n_movies}) \
      <(echo -e ~{sep="\n" qc_pass_sex}) \
      <(echo -e ~{sep="\n" qc_pass_swap}) \
      <(echo -e ~{sep="\n" qc_pass_relatedness}) \
      <(echo -e ~{sep="\n" sex}) \
      <(echo -e ~{sep="\n" coverage_mean}) \
      <(echo -e ~{sep="\n" read_count}) \
      <(echo -e ~{sep="\n" unique_read_count}) \
      <(echo -e ~{sep="\n" read_quality_mean}) \
      <(echo -e ~{sep="\n" read_quality_median}) \
      <(echo -e ~{sep="\n" read_quality_stdev}) \
      <(echo -e ~{sep="\n" read_length_mean}) \
      <(echo -e ~{sep="\n" read_length_median}) \
      <(echo -e ~{sep="\n" read_length_stdev}) \
      >> ~{cohort_id}.quality_control_summary.tsv
  >>>

  output {
    File quality_control_summary = "~{cohort_id}.quality_control_summary.tsv"
  }
  runtime {
    cpu: 1
    memory: "2 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "ubuntu:20.04"
  }
}