version 1.0

import "../structs.wdl"

task somalier_relate {

  input {
    String sample_id
    Array[File] bams
    Array[File] bam_indexes

    File reference_fasta
    File reference_index
    File sample_swap_sites_vcf

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 20
  Int threads = 4
  Int mem_gb = 4 * threads

  command<<<
    set -o pipefail

    # symlink bam and bai to same location so they can be found together
    for bam in ~{sep=" " bams}; do ln -s "$(readlink -f $bam)" . ; done
    for bai in ~{sep=" " bam_indexes}; do ln -s "$(readlink -f $bai)" . ; done

    # extract genotype-like information for each movie at selected sites 
    ls -1 ./*.bam \
    | parallel --jobs ~{threads} \
      somalier extract \
        --fasta=~{reference_fasta} \
        --sites=~{sample_swap_sites_vcf} \
        --out-dir=extracted/ \
        --sample-prefix="$(basename {} .bam)" \
        {}

    # calculate relatedness among movies from extracted, genotype-like information
    somalier relate \
      --min-depth=4 \
      --infer \
      --output-prefix=~{sample_id}.somalier \
      extracted/*.somalier
    
    # get minimum pairwise relatedness
    awk 'NR>1 {print $3}' ~{sample_id}.somalier.pairs.tsv | sort -n | head -1 > min_relatedness.txt
    # get inferred sex
    LOW=$(awk 'NR>1 {print $5}' ~{sample_id}.somalier.samples.tsv | sort -n | head -1)
    HIGH=$(awk 'NR>1 {print $5}' ~{sample_id}.somalier.samples.tsv | sort -n | tail -1)
    if [ $HIGH -eq $LOW ]; then
      if [ $HIGH -eq 1 ]; then
        echo "M" > inferred_sex.txt
      elif [ $HIGH -eq 2 ]; then
        echo "F" > inferred_sex.txt
      fi
    else
      echo "U" > inferred_sex.txt
    fi
  >>>

  output {
    File groups = "~{sample_id}.somalier.groups.tsv"
    File html = "~{sample_id}.somalier.html"
    File pairs = "~{sample_id}.somalier.pairs.tsv"
    File samples = "~{sample_id}.somalier.samples.tsv"
    Float min_relatedness = read_float("min_relatedness.txt")
    String inferred_sex = read_string("inferred_sex.txt")
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
    docker: "~{runtime_attributes.container_registry}/somalier:0.2.16"
  }
}
