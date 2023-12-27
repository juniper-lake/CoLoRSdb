version 1.0

# Calculate VCF stats

import "../structs.wdl"

task filter_norm_deepvariant {
  input {
    File vcf
    File vcf_index

    File reference_fasta

    Array[String] chromosomes

    RuntimeAttributes runtime_attributes
  }

  String vcf_basename = basename(vcf, ".vcf.gz")

  Int threads = 8
  Int disk_size = ceil((size(vcf, "GB") + size(reference_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --help

    ln -s ~{vcf}
    ln -s ~{vcf_index}

    # filter, normalize, and split multiallelics by chromosome
    parallel -j ~{threads} \
      'bcftools norm --regions {} --multiallelics - --fasta-ref ~{reference_fasta} ~{vcf_basename}.vcf.gz \
        | bcftools sort --output-type z --output {}.norm.vcf.gz' ::: ~{sep=" " chromosomes}

    # concat chromosome vcfs
    bcftools concat \
      --threads ~{threads - 1} \
      --output-type z \
      --output ~{vcf_basename}.norm.vcf.gz \
      {~{sep="," chromosomes}}.norm.vcf.gz

    bcftools index \
      --threads ~{threads - 1} \
      --tbi \
      ~{vcf_basename}.norm.vcf.gz
  >>>

  output {
    File normalized_vcf = "~{vcf_basename}.norm.vcf.gz"
    File normalized_vcf_index = "~{vcf_basename}.norm.vcf.gz.tbi"
  }

  runtime {
    cpu: threads
    memory: "16 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:1f7c3d55cac4210514a7f9108d544282fea8ad68ff54c97baa3b0a4fa2628300"
  }
}

task structural_variant_stats {
  input {
    File vcf
    Array[String] sample_ids

    Array[String] autosomes

    RuntimeAttributes runtime_attributes
  }

  String vcf_basename = basename(vcf, ".gz")

  Int threads = 2
  Int disk_size = ceil((size(vcf, "GB")) * 1.5 + 20)

  command <<<
    set -euo pipefail

    bcftools --help

    for SAMPLE in ~{sep=" " sample_ids}; do
      bcftools query \
        --targets ~{sep="," autosomes} \
        --samples $SAMPLE \
        --format '%SVTYPE\t[%GT]\n]' \
        ~{vcf} > ~{vcf_basename}.$SAMPLE.txt

      cat ~{vcf_basename}.$SAMPLE.txt \
        | datamash -s groupby 1,2 count 1 \
        | awk -v OFS="\t" -v sample=$SAMPLE '{print sample,$0}' \
        >> ~{vcf_basename}.autosomes.stats.tsv
    done
  >>>

  output {
    File stats = "~{vcf_basename}.autosomes.stats.tsv"
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
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:1f7c3d55cac4210514a7f9108d544282fea8ad68ff54c97baa3b0a4fa2628300"
  }
}

task concat_pbsv_vcfs {
  input {
    Array[File] vcfs
    String output_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads = 4
  Int mem_gb = 4
  Int disk_size = ceil(size(vcfs, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    if [[ "~{length(vcfs)}" -eq 1 ]]; then
      cp ~{vcfs[0]} concat.tmp.vcf
    else
      # will throw error if vcfs are not in correct order
      # could --allow-overlaps but would require sort after and that's compute intensive
      bcftools concat \
        --output concat.tmp.vcf \
        ~{sep=' ' vcfs}
    fi

    # keep pass-only variants
    bcftools view \
      --output-type z \
      --apply-filters PASS \
      --output concat.tmp.vcf.gz \
      concat.tmp.vcf

    # must change AD values to DV/DR for jasmine compatibility
    bcftools query \
      --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD{1}][\t%AD{0}]\n' \
      concat.tmp.vcf.gz \
      | bgzip --stdout \
      > annot.txt.gz

    tabix \
      --sequence 1 \
      --begin 2 \
      --end 2 \
      annot.txt.gz

    echo -e '##FORMAT=<ID=DV,Number=1,Type=String,Description="The number of reads supporting the variant sequence">' \
      > hdr.txt
    echo -e '##FORMAT=<ID=DR,Number=1,Type=String,Description="The number of reads supporting the reference sequence">' \
      >> hdr.txt

    bcftools annotate \
      --annotations annot.txt.gz \
      --header-lines hdr.txt \
      --columns CHROM,POS,REF,ALT,FORMAT/DV,FORMAT/DR \
      --output-type v \
      --output ~{output_prefix}.vcf \
      concat.tmp.vcf.gz

    bgzip \
      --threads ~{threads} \
      --stdout \
      ~{output_prefix}.vcf \
      > ~{output_prefix}.vcf.gz

    bcftools index \
      --threads ~{threads - 1} \
      --tbi \
      ~{output_prefix}.vcf.gz
  >>>

  output {
    File concatenated_zipped_vcf = "~{output_prefix}.vcf.gz"
    File concatenated_zipped_vcf_index = "~{output_prefix}.vcf.gz.tbi"
    File concatenated_vcf = "~{output_prefix}.vcf"
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
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:1f7c3d55cac4210514a7f9108d544282fea8ad68ff54c97baa3b0a4fa2628300"
  }
}

task merge_vcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    String output_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads = 4
  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

command {
    set -euo pipefail

    # increase open file limit
    ulimit -Sn 65536

    bcftools merge \
      --merge none \
      --threads ~{threads - 1} \
      --output-type z \
      ~{sep=" " vcfs} \
      > ~{output_prefix}.vcf.gz

    bcftools index \
      --threads ~{threads - 1} \
      --tbi \
      ~{output_prefix}.vcf.gz
  }

  output {
    File merged_vcf = "~{output_prefix}.vcf.gz"
    File merged_vcf_index = "~{output_prefix}.vcf.gz.tbi"
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
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:1f7c3d55cac4210514a7f9108d544282fea8ad68ff54c97baa3b0a4fa2628300"
  }
}

task filter_zip_index_sniffles {
  input {
    File vcf
    Array[String] chromosomes
    RuntimeAttributes runtime_attributes
  }

  String vcf_basename = basename(vcf)
  Int threads = 4
  Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --help

    bcftools view \
      --threads ~{threads - 1} \
      --targets ~{sep="," chromosomes} \
      --apply-filters PASS \
      --output-type z \
      --output ~{vcf_basename}.gz \
      ~{vcf}

    bcftools index \
      --threads ~{threads - 1} \
      --tbi \
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
    disks: "local-disk ~{disk_size} LOCAL"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:1f7c3d55cac4210514a7f9108d544282fea8ad68ff54c97baa3b0a4fa2628300"
  }
}

task reheader_zip_index_jasminesv {
  input {
    File vcf
    Array[String] sample_ids
    RuntimeAttributes runtime_attributes
  }

  String vcf_basename = basename(vcf)
  Int threads = 4
  Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    echo -e "~{sep="\n" sample_ids}" > samples.txt

    bcftools --help

    bcftools reheader \
      --threads ~{threads - 1} \
      --samples samples.txt \
      ~{vcf} \
    | bcftools sort \
      --output-type z \
      --output ~{vcf_basename}.gz

    bcftools index \
      --threads ~{threads - 1} \
      --tbi \
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
    disks: "local-disk ~{disk_size} LOCAL"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:1f7c3d55cac4210514a7f9108d544282fea8ad68ff54c97baa3b0a4fa2628300"
  }
}
