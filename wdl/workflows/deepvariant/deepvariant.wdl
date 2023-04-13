version 1.0

import "../../structs.wdl"

workflow deepvariant {

  input {
    String sample_id
    Array[IndexData] aligned_bams

    String reference_name
    IndexData reference_fasta

    String deepvariant_version

    RuntimeAttributes default_runtime_attributes
  }

  Int shards = 64

  scatter (bam_object in aligned_bams) {
		File bam = bam_object.data
		File bam_index = bam_object.index
	}

  call make_examples {
    input: 
      sample_id = sample_id,
      bams = bam,
      bam_indexes = bam_index,
      reference_fasta = reference_fasta.data,
      reference_index = reference_fasta.index,
      threads = shards,
      deepvariant_version = deepvariant_version,
      runtime_attributes = default_runtime_attributes
  }
  
  call call_variants {
    input: 
      sample_id = sample_id,
      reference_name = reference_name,
      example_tfrecords = make_examples.example_tfrecords,
      threads = shards,
      deepvariant_version = deepvariant_version,
      runtime_attributes = default_runtime_attributes
  }

  call postprocess_variants {
    input:
      sample_id = sample_id,
      tfrecord = call_variants.tfrecord,
      nonvariant_site_tfrecords = make_examples.nonvariant_site_tfrecords,
      reference_name = reference_name,
      reference_fasta = reference_fasta.data,
      reference_index = reference_fasta.index,
      shards = shards,
      deepvariant_version = deepvariant_version,
      runtime_attributes = default_runtime_attributes

  }

  output {
    IndexData vcf = {"data": postprocess_variants.vcf, "index": postprocess_variants.vcf_index}
    IndexData gvcf = {"data": postprocess_variants.gvcf, "index": postprocess_variants.gvcf_index}
    File report = postprocess_variants.report
  }


  
}


task make_examples {

  input {
    String sample_id
    Array[File] bams
    Array[File] bam_indexes

    File reference_fasta
    File reference_index

    Int threads
    String deepvariant_version

    RuntimeAttributes runtime_attributes
  }
  
  Int mem_gb = 4 * threads
	Int disk_size = ceil(size(bams[0], "GB") * length(bams) * 2 + 50)

  command <<<
    set -euo pipefail
    
    seq 0 ~{threads - 1} \
    | parallel --jobs ~{threads} \
      /opt/deepvariant/bin/make_examples \
        --norealign_reads \
        --vsc_min_fraction_indels 0.12 \
        --pileup_image_width 199 \
        --track_ref_reads \
        --phase_reads \
        --partition_size=25000 \
        --max_reads_per_partition=600 \
        --alt_aligned_pileup=diff_channels \
        --add_hp_channel \
        --sort_by_haplotypes \
        --parse_sam_aux_fields \
        --min_mapping_quality=1 \
        --mode calling \
        --ref ~{reference_fasta} \
        --reads ~{sep="," bams} \
        --examples ~{sample_id}.examples.tfrecord@~{threads}.gz \
        --gvcf ~{sample_id}.gvcf.tfrecord@~{threads}.gz \
        --task {}
  >>>

  output {
    Array[File] example_tfrecords = glob("~{sample_id}.examples.tfrecord*.gz")
    Array[File] nonvariant_site_tfrecords = glob("~{sample_id}.gvcf.tfrecord*.gz")
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
    docker: "google/deepvariant:~{deepvariant_version}"
  }
}


task call_variants {
  
  input {
    String sample_id
    String reference_name
    Array[File] example_tfrecords

    String deepvariant_version
    Int threads

    RuntimeAttributes runtime_attributes
  }
  
  String example_tfrecord_path = sub(example_tfrecords[0], "/" + basename(example_tfrecords[0]), "")
  String outfile = "~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz"
	Int disk_size = ceil(size(example_tfrecords[0], "GB") * length(example_tfrecords) * 2 + 100)

  command {
    set -euo pipefail
    
    /opt/deepvariant/bin/call_variants \
      --outfile ~{outfile} \
      --examples ~{example_tfrecord_path}/~{sample_id}.examples.tfrecord@~{threads}.gz \
      --checkpoint "/opt/models/pacbio/model.ckpt"
  }

  output {
    File tfrecord = outfile
  }

  runtime {
    cpu: threads
    memory: "8 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "google/deepvariant:~{deepvariant_version}"
  }
}


task postprocess_variants {
  
  input {
    String sample_id
    File tfrecord
    Array[File] nonvariant_site_tfrecords
    
    String reference_name
    File reference_fasta
    File reference_index

    String deepvariant_version
    Int shards

    RuntimeAttributes runtime_attributes
  }

  String nonvariant_site_tfrecord_path = sub(nonvariant_site_tfrecords[0], "/" + basename(nonvariant_site_tfrecords[0]), "")
  Int disk_size = ceil((size(tfrecord, "GB") + size(reference_fasta, "GB") + size(nonvariant_site_tfrecords[0], "GB") * length(nonvariant_site_tfrecords)) * 2 + 20)
  String outfile = "~{sample_id}.~{reference_name}.deepvariant.vcf.gz"
  String gvcf_outfile = "~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz"

  command {
    set -euo pipefail
    
    /opt/deepvariant/bin/postprocess_variants \
      --ref ~{reference_fasta} \
      --infile ~{tfrecord} \
      --outfile ~{outfile} \
      --nonvariant_site_tfrecord_path ~{nonvariant_site_tfrecord_path}/~{sample_id}.gvcf.tfrecord@~{shards}.gz \
      --gvcf_outfile ~{gvcf_outfile}
  }

  output {
    File vcf = outfile
    File vcf_index = "~{outfile}.tbi"
    File gvcf = gvcf_outfile
    File gvcf_index = "~{gvcf_outfile}.tbi"
    File report = "~{sample_id}.~{reference_name}.deepvariant.visual_report.html"
  } 

  runtime {
    cpu: 1
    memory: "32 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "google/deepvariant:~{deepvariant_version}"
  }
}
