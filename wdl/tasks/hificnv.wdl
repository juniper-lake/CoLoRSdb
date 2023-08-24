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
      --output-prefix hificnv

    # because hificnv outputs symbolic alleles only (<DEL> and <DUP>), CNVs starting at the same position
    # but possessing different length get merged and the QUAL + INFO fields are lost
    # therefore we move QUAL to FORMAT/CNQ and add an ID field to each variant on which to merge

    bcftools --version
    bgzip --version

    # extract QUAL and make ID fields from VCF
    bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\thificnv.%INFO/SVTYPE\.%CHROM\.%POS\.%INFO/END\n' \
      hificnv.~{sample_id}.vcf.gz \
      | bgzip -c \
      > annot.txt.gz

    tabix --version

    # index the file with tabix
    tabix -s1 -b2 -e2 annot.txt.gz

    # create new header line for QUAL annotation
    echo -e '##FORMAT=<ID=CNQ,Number=1,Type=Integer,Description="Average of the next-most-likely copy-number state for each bin from the HMM">' >> hdr.txt

    # annotate VCF with new info
    bcftools annotate \
      -x QUAL \
      -a annot.txt.gz \
      -h hdr.txt \
      -c CHROM,POS,REF,ALT,FORMAT/CNQ,ID \
      -o ~{sample_id}.~{reference_name}.hificnv.vcf \
      hificnv.~{sample_id}.vcf.gz

    # compress 
    bgzip \
      --threads ~{threads} \
      ~{sample_id}.~{reference_name}.hificnv.vcf
    
    tabix \
      --preset vcf \
      ~{sample_id}.~{reference_name}.hificnv.vcf.gz

  >>>

  output {
    File cnv_vcf = "~{sample_id}.~{reference_name}.hificnv.vcf.gz"
    File cnv_vcf_index = "~{sample_id}.~{reference_name}.hificnv.vcf.gz.tbi"
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
    docker: "~{runtime_attributes.container_registry}/hificnv@sha256:c8b164686c7e21cf4ecdc3b56f3c3303471df0582e409682bd824bd7a99db96a"
  }
}


task merge_hificnv_vcfs {
  input {
    Array[File] cnv_vcfs
    Array[File] cnv_vcf_indexes
    String cohort_id
    String reference_name

    RuntimeAttributes runtime_attributes
  }
  
  Int mem_gb = 8
  Int threads = 4
  Int disk_size = ceil((size(cnv_vcfs, "GB")) * 2 + 20)

  command {
    set -euo pipefail

    # increase open file limit
    ulimit -Sn 65536
    
    bcftools --version

    bcftools merge \
      --merge id \
      ~{sep=" " cnv_vcfs} \
      > ~{cohort_id}.~{reference_name}.hificnv.vcf
    
    bgzip --version

    bgzip \
      --threads ~{threads} \
      ~{cohort_id}.~{reference_name}.hificnv.vcf
    
    tabix \
      --preset vcf \
      ~{cohort_id}.~{reference_name}.hificnv.vcf.gz
  }

  output {
    File merged_cnv_vcf = "~{cohort_id}.~{reference_name}.hificnv.vcf.gz"
    File merged_cnv_vcf_index = "~{cohort_id}.~{reference_name}.hificnv.vcf.gz.tbi"
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
    docker: "~{runtime_attributes.container_registry}/hificnv@sha256:c8b164686c7e21cf4ecdc3b56f3c3303471df0582e409682bd824bd7a99db96a"
  }
}
