version 1.0

# Quality control of samples (sex, sample swaps, related samples) using somalier

import "../structs.wdl"

task somalier_sample_swap {
  input {
    String sample_id
    Array[File] bams
    Array[File] bam_indexes
    Array[String] movie_names

    File reference_fasta
    File reference_index
    File somalier_sites_vcf

    RuntimeAttributes runtime_attributes
  }

  Int n_files = length(bams)
  Int disk_size = ceil((size(bams, "GB") + size(reference_fasta, "GB") + size(somalier_sites_vcf, "GB")) * 2.5 + 5)
  Int threads = n_files

  command<<<
    set -euo pipefail

    export SOMALIER_REPORT_ALL_PAIRS=1

    # get minimum pairwise relatedness (output 1 if there's only one movie)
    if [ ~{n_files} -eq 1 ]; then
      echo "1.0" > min_relatedness.txt
      echo "0" > min_relatedness_n_sites.txt
    else
      paste <(echo -e "~{sep="\n" bams}") <(echo -e "~{sep="\n" movie_names}") > parallel_variables.txt
      mkdir extracted
      parallel \
        --jobs ~{threads} \
        --colsep '\t' \
        --arg-file parallel_variables.txt \
        'somalier extract \
          --fasta=~{reference_fasta} \
          --sites=~{somalier_sites_vcf} \
          --out-dir={2} \
          --sample-prefix={2} \
          {1} && mv {2}/~{sample_id}.somalier extracted/{2}.somalier'

      # calculate relatedness among movies
      somalier relate \
        --min-depth=5 \
        --min-ab=0.1 \
        --output-prefix=~{sample_id}.somalier \
        extracted/*.somalier

      # get line corresponding to pairwise movie comparison with minimum relatedness
      LINE=$(grep -v "^#" ~{sample_id}.somalier.pairs.tsv | sort -k3 -g | head -n1)
      echo $LINE | cut -f3 -d " " > min_relatedness.txt
      echo $LINE | cut -f14 -d " " > min_relatedness_n_sites.txt
    fi
  >>>

  output {
    File? groups = "~{sample_id}.somalier.groups.tsv"
    File? html = "~{sample_id}.somalier.html"
    File? pairs = "~{sample_id}.somalier.pairs.tsv"
    File? samples = "~{sample_id}.somalier.samples.tsv"
    Float min_relatedness = read_float("min_relatedness.txt")
    Int min_relatedness_n_sites = read_int("min_relatedness_n_sites.txt")
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
    docker: "~{runtime_attributes.container_registry}/somalier@sha256:91fa892557e19a60de67e99f0a192aa2eccfec8eb2a1af05427ae48d81a24cc0"
  }
}

task somalier_extract {
  input {
    String sample_id
    String? sample_prefix
    File bam
    File bam_index

    File reference_fasta
    File reference_index
    File somalier_sites_vcf

    RuntimeAttributes runtime_attributes
  }

  String out_dir = select_first([sample_prefix, "extracted"])
  Int disk_size = ceil((size(bam, "GB") + size(reference_fasta, "GB")) * 2.5 + 5)
  Int threads = 1

  command<<<
    set -euo pipefail

    # extract sites
    somalier extract \
      --fasta=~{reference_fasta} \
      --sites=~{somalier_sites_vcf} \
      --out-dir=~{out_dir} \
      ~{"--sample-prefix=" + sample_prefix} \
      ~{bam}
  >>>

  output {
    File extracted_sites = "~{out_dir}/~{sample_id}.somalier"
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
    docker: "~{runtime_attributes.container_registry}/somalier@sha256:91fa892557e19a60de67e99f0a192aa2eccfec8eb2a1af05427ae48d81a24cc0"
  }
}

task somalier_relate_samples {
  input {
    String cohort_id
    String reference_name
    Array[File] extracted_sites
    Array[String] sample_ids
    Array[Float] coverages
    Float max_pairwise_relatedness

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 20 # might want to increase this
  Int threads = 1

  command<<<
    set -euo pipefail

    export SOMALIER_REPORT_ALL_PAIRS=1

    # increase open file limit
    ulimit -Sn 65536

    # calculate relatedness among samples
    somalier relate \
      --min-depth=6 \
      --min-ab=0.3 \
      --output-prefix=~{cohort_id}.~{reference_name}.somalier \
      --infer \
      ~{sep=" " extracted_sites}

    # find samples that have relatedness > max_pairwise_relatedness
    screen_related_samples.py \
      ~{cohort_id}.~{reference_name}.somalier.pairs.tsv \
      --max_relatedness ~{max_pairwise_relatedness} \
      --sample_order ~{sep=" " sample_ids} \
      --coverages ~{sep=" " coverages} \
      --outfile ~{cohort_id}.related_samples_to_remove.tsv

    for SAMPLE_ID in ~{sep=" " sample_ids}; do
      awk -v sample=$SAMPLE_ID '$1==sample {print $2}' ~{cohort_id}.related_samples_to_remove.tsv >> keep_drop.txt
      awk -v sample=$SAMPLE_ID '$1==sample {print $3}' ~{cohort_id}.related_samples_to_remove.tsv >> n_relations.txt
      SEX=$(awk -v sample=$SAMPLE_ID '$2==sample {print $5}' ~{cohort_id}.~{reference_name}.somalier.samples.tsv)
      Y_N=$(awk -v sample=$SAMPLE_ID '$2==sample {print $25}' ~{cohort_id}.~{reference_name}.somalier.samples.tsv)
      X_N=$(awk -v sample=$SAMPLE_ID '$2==sample {print $20}' ~{cohort_id}.~{reference_name}.somalier.samples.tsv)
      X_HET=$(awk -v sample=$SAMPLE_ID '$2==sample {print $22}' ~{cohort_id}.~{reference_name}.somalier.samples.tsv)
      X_HET_PROP="$((100*X_HET/X_N))"
      if [ $SEX -eq 1 ]; then
        echo "male" >> inferred_sex.txt
      elif [ $SEX -eq 2 ]; then
        echo "female" >> inferred_sex.txt
      else
        if [ $Y_N -gt 0 ] && [ $X_HET_PROP -lt 2 ]; then
          echo "male" >> inferred_sex.txt
        elif [ $Y_N -lt 1 ] && [ $X_HET_PROP -gt 10 ]; then
          echo "female" >> inferred_sex.txt
        else
          echo "unknown" >> inferred_sex.txt
        fi
      fi
    done
  >>>

  output {
    File? groups = "~{cohort_id}.~{reference_name}.somalier.groups.tsv"
    File html = "~{cohort_id}.~{reference_name}.somalier.html"
    File pairs = "~{cohort_id}.~{reference_name}.somalier.pairs.tsv"
    File samples = "~{cohort_id}.~{reference_name}.somalier.samples.tsv"
    Array[String] qc_related_keep_drop = read_lines("keep_drop.txt")
    Array[String] n_relations = read_lines("n_relations.txt")
    Array[String] inferred_sexes = read_lines("inferred_sex.txt")
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
    docker: "~{runtime_attributes.container_registry}/somalier@sha256:91fa892557e19a60de67e99f0a192aa2eccfec8eb2a1af05427ae48d81a24cc0"
  }
}
