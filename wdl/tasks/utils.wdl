version 1.0

# Read sample ids and movie paths from file

import "../structs.wdl"

task unzip_reference_bundle {
  input {
    File reference_bundle

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(reference_bundle) * 5) + 10

  command <<<
    set -euo pipefail

    tar -xzf ~{reference_bundle}
  >>>

  output {
    ReferenceData grch38 = object {
      "name": "GRCh38",
      "fasta": "colorsdb_resources/GRCh38/human_GRCh38_no_alt_analysis_set.fasta",
      "fasta_index": "colorsdb_resources/GRCh38/human_GRCh38_no_alt_analysis_set.fasta.fai",
      "chromosomes": [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
        "chrM"
    ],
    "autosomes": [
      "chr1",
      "chr2",
      "chr3",
      "chr4",
      "chr5",
      "chr6",
      "chr7",
      "chr8",
      "chr9",
      "chr10",
      "chr11",
      "chr12",
      "chr13",
      "chr14",
      "chr15",
      "chr16",
      "chr17",
      "chr18",
      "chr19",
      "chr20",
      "chr21",
      "chr22"
    ],
    "non_diploid_regions": "colorsdb_resources/GRCh38/vcfparser.GRCh38.ploidy.txt",
    "hificnv_exclude_bed": "colorsdb_resources/GRCh38/hificnv.cnv.excluded_regions.hg38.bed.gz",
    "hificnv_exclude_bed_index": "colorsdb_resources/GRCh38/hificnv.cnv.excluded_regions.hg38.bed.gz.tbi",
    "hificnv_expected_bed_male": "colorsdb_resources/GRCh38/hificnv.male_expected_cn.hg38.bed",
    "hificnv_expected_bed_female": "colorsdb_resources/GRCh38/hificnv.female_expected_cn.hg38.bed",
    "tandem_repeat_bed": "colorsdb_resources/GRCh38/human_GRCh38_no_alt_analysis_set.trf.bed",
    "trgt_tandem_repeat_beds": [
      "colorsdb_resources/GRCh38/trgt.adotto_repeats.hg38.bed",
      "colorsdb_resources/GRCh38/trgt.pathogenic_repeats.hg38.bed",
      "colorsdb_resources/GRCh38/trgt.repeat_catalog.hg38.bed"
    ],
    "somalier_sites_vcf": "colorsdb_resources/GRCh38/somalier.sites.hg38.vcf.gz",
    "peddy_sites": "colorsdb_resources/GRCh38/peddy.GRCH38.sites",
    "peddy_bin": "colorsdb_resources/GRCh38/peddy.GRCH38.sites.bin.gz"
    }

  ReferenceData chm13 = object {
    "name": "CHM13",
    "fasta": "colorsdb_resources/CHM13/human_chm13v2.0_maskedY_rCRS.fasta",
    "fasta_index": "colorsdb_resources/CHM13/human_chm13v2.0_maskedY_rCRS.fasta.fai",
    "chromosomes": [
      "chr1",
      "chr2",
      "chr3",
      "chr4",
      "chr5",
      "chr6",
      "chr7",
      "chr8",
      "chr9",
      "chr10",
      "chr11",
      "chr12",
      "chr13",
      "chr14",
      "chr15",
      "chr16",
      "chr17",
      "chr18",
      "chr19",
      "chr20",
      "chr21",
      "chr22",
      "chrX",
      "chrY",
      "chrM"
    ],
    "autosomes": [
      "chr1",
      "chr2",
      "chr3",
      "chr4",
      "chr5",
      "chr6",
      "chr7",
      "chr8",
      "chr9",
      "chr10",
      "chr11",
      "chr12",
      "chr13",
      "chr14",
      "chr15",
      "chr16",
      "chr17",
      "chr18",
      "chr19",
      "chr20",
      "chr21",
      "chr22"
    ],
    "non_diploid_regions": "colorsdb_resources/CHM13/vcfparser.CHM13.ploidy.txt",
    "tandem_repeat_bed": "colorsdb_resources/CHM13/human_chm13v2.0_maskedY_rCRS.trf.bed",
    "somalier_sites_vcf": "colorsdb_resources/CHM13/somalier.sites.chm13v2.T2T.vcf.gz"
    }
  }

  runtime {
    cpu: 1
    memory: "1 GB"
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

task read_sample_sheet {
  input {
    File sample_sheet

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(sample_sheet)) + 5

  command <<<
    set -euo pipefail

    grep -v "^#" ~{sample_sheet} | grep -v -e '^$' > temp_sample_sheet.tsv
    cut -f1 temp_sample_sheet.tsv > sample_ids.txt
    cut -f2 temp_sample_sheet.tsv | sed 's/,/\t/g' > movies.tsv
    cut -f3 temp_sample_sheet.tsv > qc_pass.txt
    cut -f4 temp_sample_sheet.tsv > sexes.txt

    # are the number of rows the same for each column
    if [ "$(grep -v -e '^$' sample_ids.txt | wc -l)" -ne "$(grep -v -e '^$' movies.tsv | wc -l)" ]; then
      echo "The number of sample id column and movies column in the sample sheet do not have the same number of rows." 1>&2
      exit 1
    fi

    # are all sample ids unique
    if [ "$(sort sample_ids.txt | uniq -d)" ]; then
      echo "The sample ids in the sample sheet are not unique." 1>&2
      exit 1
    fi

    # are all movie paths unique
    if [ "$(cat movies.tsv | sed 's/\t\t*/\n/g' | sort | uniq -d)" ]; then
      echo "The movie paths in the sample sheet are not unique." 1>&2
      exit 1
    fi
  >>>

  output {
    Array[String] sample_ids = read_lines("sample_ids.txt")
    Array[Array[String]] movies = read_tsv("movies.tsv")
    Array[String] qc_pass = read_lines("qc_pass.txt")
    Array[String] sexes = read_lines("sexes.txt")
  }

  runtime {
    cpu: 1
    memory: "1 GB"
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

task summarize_qc {
  input {
    String cohort_id
    String reference_name
    Array[String] sample_ids
    Array[Float] movie_relatedness_qc_threshold
    Array[Float] sample_relatedness_qc_threshold
    Array[Boolean] qc_pass_swap
    Array[Boolean] qc_pass_sex
    Array[Boolean] qc_pass_relatedness
    Array[Boolean] qc_pass_bams
    Array[Boolean] qc_pass_combined
    Array[Float] min_movie_relatedness
    Array[Int] min_movie_relatedness_n_sites
    Array[String] n_relations
    Array[String] sex
    Array[Int] n_movies
    Array[Float] coverage_mean
    Array[Int] read_count
    Array[Int] unique_read_count
    Array[Float] read_quality_mean
    Array[Int] read_quality_median
    Array[Float] read_quality_stdev
    Array[Float] read_length_mean
    Array[Int] read_length_median
    Array[Float] read_length_stdev

    Array[Array[String]] peek_a_bam_tsv
    Array[String] peek_a_bam_fields = ["sample_id","file","samples","potentially_multiplexed","demultiplexed","ccs","hifi","kinetics","base_modification","aligned","haplotagged"]

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 20

  command <<<
    set -euo pipefail

    paste <(echo -e "sample_id\n~{sep="\n" sample_ids}") \
      <(echo -e "movie_relatedness_qc_threshold\n~{sep="\n" movie_relatedness_qc_threshold}") \
      <(echo -e "sample_relatedness_qc_threshold\n~{sep="\n" sample_relatedness_qc_threshold}") \
      <(echo -e "qc_pass_swap\n~{sep="\n" qc_pass_swap}") \
      <(echo -e "qc_pass_relatedness\n~{sep="\n" qc_pass_relatedness}") \
      <(echo -e "qc_pass_bams\n~{sep="\n" qc_pass_bams}") \
      <(echo -e "qc_pass_sex\n~{sep="\n" qc_pass_sex}") \
      <(echo -e "qc_pass_combined\n~{sep="\n" qc_pass_combined}") \
      <(echo -e "min_movie_relatedness\n~{sep="\n" min_movie_relatedness}") \
      <(echo -e "min_movie_relatedness_n_sites\n~{sep="\n" min_movie_relatedness_n_sites}") \
      <(echo -e "n_relations\n~{sep="\n" n_relations}") \
      <(echo -e "sex\n~{sep="\n" sex}") \
      <(echo -e "n_movies\n~{sep="\n" n_movies}") \
      <(echo -e "coverage_mean\n~{sep="\n" coverage_mean}") \
      <(echo -e "read_count\n~{sep="\n" read_count}") \
      <(echo -e "unique_read_count\n~{sep="\n" unique_read_count}") \
      <(echo -e "read_quality_mean\n~{sep="\n" read_quality_mean}") \
      <(echo -e "read_quality_median\n~{sep="\n" read_quality_median}") \
      <(echo -e "read_quality_stdev\n~{sep="\n" read_quality_stdev}") \
      <(echo -e "read_length_mean\n~{sep="\n" read_length_mean}") \
      <(echo -e "read_length_median\n~{sep="\n" read_length_median}") \
      <(echo -e "read_length_stdev\n~{sep="\n" read_length_stdev}") \
      > ~{cohort_id}.~{reference_name}.quality_control_summary.tsv

    cat <(echo -e "~{sep="\t" peek_a_bam_fields}") ~{write_tsv(peek_a_bam_tsv)} \
      > ~{cohort_id}.~{reference_name}.quality_control_bams.tsv
  >>>

  output {
    File sample_qc_tsv = "~{cohort_id}.~{reference_name}.quality_control_summary.tsv"
    File bam_qc_tsv = "~{cohort_id}.~{reference_name}.quality_control_bams.tsv"

  }
  runtime {
    cpu: 1
    memory: "1 GB"
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
