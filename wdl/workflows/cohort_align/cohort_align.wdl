version 1.0

import "../../tasks/pbmm2.wdl" as Pbmm2
import "../../tasks/utils.wdl" as Utils
import "../../tasks/mosdepth.wdl" as Mosdepth

workflow cohort_align {

  input {
    String cohort_id
    Array[Sample] samples

    ReferenceData reference

    RuntimeAttributes default_runtime_attributes
  }

  scatter (sample in samples) {
    # align each movie with pbmm2
    scatter (idx in range(length(sample.movies))) {    

      # get movie name from movie path but append "_movie1", "_movie2", etc in case file names
      # are not unique or don't follow HiFi naming convention
      String movie_name = sub(basename(sample.movies[idx]), "\\..*", "_movie~{idx}")
      call Pbmm2.pbmm2 {
          input: 
            movie = sample.movies[idx],
            out_prefix = "~{sample.sample_id}.~{movie_name}",
            sample_id = sample.sample_id,
            reference_name = reference.name,
            reference_fasta = reference.fasta.data,
            reference_index = reference.fasta.index,
            runtime_attributes = default_runtime_attributes
      }
    }

    # merge aligned bams
    call Utils.merge_bams {
      input:
        bams = pbmm2.aligned_bam,
        output_bam_name = "~{sample.sample_id}.~{reference.name}.bam",
        runtime_attributes = default_runtime_attributes
    }

    IndexData merged_bam = {
      "data": merge_bams.merged_bam, 
      "index": merge_bams.merged_bam_index
    }

    # calculate depth of merged bam
    call Mosdepth.mosdepth {
      input:
        aligned_bam = merged_bam.data,
        aligned_bam_index = merged_bam.index,
        runtime_attributes = default_runtime_attributes
    }
    
    String sample_id = sample.sample_id
    Int n_movies = length(sample.movies)
    
    AlignedSample aligned_sample = object {
      "sample_id": sample.sample_id,
      "aligned_bam": merged_bam,
      "sex": "UNKNOWN",
      "qc_pass": false
    }
  }

  call summarize_coverage {
    input:
      cohort_id = cohort_id,
      reference_name = reference.name,
      sample_ids = sample_id,
      n_movies = n_movies,
      coverage_mean = mosdepth.mean_coverage,
      runtime_attributes = default_runtime_attributes
  }

  output {
    Array[AlignedSample] aligned_samples = aligned_sample
    File coverage_summary_tsv = summarize_coverage.coverage_summary
  }
}


task summarize_coverage {
  input {
    String cohort_id
    String reference_name
    Array[String] sample_ids
    Array[Int] n_movies
    Array[Float] coverage_mean
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 20

  command <<<
    set -euo pipefail
    
    paste <(echo -e "sample_id\n~{sep="\n" sample_ids}") \
      <(echo -e "n_movies\n~{sep="\n" n_movies}") \
      <(echo -e "coverage_mean\n~{sep="\n" coverage_mean}") \
      > ~{cohort_id}.quality_control_summary.tsv
  >>>

  output {
    File coverage_summary = "~{cohort_id}.~{reference_name}.coverage_summary.tsv"
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