version 1.0

# Extract genotypes from known loci and compare across samples/movies for QC

import "../structs.wdl"

task somalier_extract {

  input {
    String sample_id
    String sample_prefix
    File bam
    File bam_index

    File reference_fasta
    File reference_index
    File somalier_sites_vcf

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(bam, "GB") + size(reference_fasta, "GB")) * 2.5 + 5)
  Int threads = 1

  command<<<

    set -euo pipefail

    # extract sites with movie prefix to check sample swaps between movies
    somalier extract \
      --fasta=~{reference_fasta} \
      --sites=~{somalier_sites_vcf} \
      --out-dir=extracted_movies \
      --sample-prefix=~{sample_prefix} \
      ~{bam}
    
    # extract sites again without movie prefix to check sample relatedness
    somalier extract \
      --fasta=~{reference_fasta} \
      --sites=~{somalier_sites_vcf} \
      --out-dir=extracted_samples \
      ~{bam}

  >>>

  output {
    File extracted_sites_per_movie = "extracted_movies/~{sample_id}.somalier"
    File extracted_sites_per_sample = "extracted_samples/~{sample_id}.somalier"
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
    docker: "~{runtime_attributes.container_registry}/somalier:0.2.16"
  }
}

task somalier_relate_movies {

  input {
    String sample_id
    Array[File] extracted_somalier_sites

    RuntimeAttributes runtime_attributes
  }

  Int n_files = length(extracted_somalier_sites)
  Int disk_size = 20
  Int threads = 1

  command<<<
    set -euo pipefail

    # calculate relatedness among movies from extracted, genotype-like information
    somalier relate \
      --min-depth=4 \
      --infer \
      --output-prefix=~{sample_id}.somalier \
      ~{sep=" " extracted_somalier_sites}
    
    # get minimum pairwise relatedness (output 1 if there's only one movie)
    if [ ~{n_files} -eq 1 ]; then
      echo "1.0" > min_relatedness.txt
    else
      awk 'NR>1 {print $3}' ~{sample_id}.somalier.pairs.tsv | sort -n | head -1 > min_relatedness.txt
    fi

    # get inferred sex
    LOW=$(awk 'NR>1 {print $5}' ~{sample_id}.somalier.samples.tsv | sort -n | head -1)
    HIGH=$(awk 'NR>1 {print $5}' ~{sample_id}.somalier.samples.tsv | sort -n | tail -1)
    if [ $HIGH -eq $LOW ]; then
      if [ $HIGH -eq 1 ]; then
        echo "MALE" > inferred_sex.txt
      elif [ $HIGH -eq 2 ]; then
        echo "FEMALE" > inferred_sex.txt
      fi
    else
      echo "UNKNOWN" > inferred_sex.txt
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
    memory: "8 GB"
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


task somalier_relate_samples {

  input {
    String cohort_id
    Array[File] extracted_somalier_sites
    Array[String] sample_ids
    Array[Float] coverages
    Float max_pairwise_relatedness

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 20 # might want to increase this
  Int threads = 1

  command<<<
    set -euo pipefail

    # calculate relatedness among movies from extracted, genotype-like information
    somalier relate \
      --min-depth=4 \
      --output-prefix=~{cohort_id}.somalier \
      ~{sep=" " extracted_somalier_sites}

    # find samples that have relatedness > 0.125
    screen_related_samples.py \
      ~{cohort_id}.somalier.pairs.tsv \
      --max_relatedness ~{max_pairwise_relatedness} \
      --sample_order ~{sep=" " sample_ids} \
      --coverages ~{sep=" " coverages} \
      --outfile ~{cohort_id}.related_samples_to_remove.tsv

    
  >>>

  output {
    File groups = "~{cohort_id}.somalier.groups.tsv"
    File html = "~{cohort_id}.somalier.html"
    File pairs = "~{cohort_id}.somalier.pairs.tsv"
    File samples = "~{cohort_id}.somalier.samples.tsv"
    Array[Array[String]] qc_related_keep_drop = read_tsv("~{cohort_id}.related_samples_to_remove.tsv")
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
    docker: "~{runtime_attributes.container_registry}/somalier:0.2.16"
  }
}
