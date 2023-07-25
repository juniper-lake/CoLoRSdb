version 1.0

# Estimate ancestry with peddy

import "../structs.wdl"

task peddy {
  input {
    String cohort_id
    Array[String] sample_ids

    File vcf
    File vcf_index

    File peddy_sites
    File peddy_bin

    RuntimeAttributes runtime_attributes
  }

  Int threads = 4
  Int mem_gb = ceil(threads * 4)
  Int disk_size = ceil((size(vcf, "GB") + size(peddy_sites, "GB") + size(peddy_bin, "GB")) * 2.5 + 10)

  command <<<
    set -euo pipefail

    # make ped file
    for SAMPLE_ID in ~{sep=" " sample_ids}; do
      echo -e "${SAMPLE_ID}\t${SAMPLE_ID}\t0\t0\t0\t0" >> ~{cohort_id}.ped
    done
    
    # enforce naming requirement for peddy to run pca 
    # binfile = sitesfile+".bin.gz"
    ln -s ~{peddy_sites} peddy.sites
    ln -s ~{peddy_bin} peddy.sites.bin.gz

    peddy --version

    peddy \
      --procs ~{threads} \
      --sites peddy.sites \
      --prefix ~{cohort_id} \
      ~{vcf} \
      ~{cohort_id}.ped
  >>>

  output {
    File het_check = "~{cohort_id}.het_check.csv"
    File sex_check = "~{cohort_id}.sex_check.csv"
    File ped_check = "~{cohort_id}.ped_check.csv"
    File background_pca = "~{cohort_id}.background_pca.json"
    File html = "~{cohort_id}.html"
    File ped = "~{cohort_id}.peddy.ped"
    File vs_html = "~{cohort_id}.vs.html"
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
    docker: "~{runtime_attributes.container_registry}/peddy@sha256:a0224eb1161a87dad49b671d5e4126fcceac3c33b1f9dffa9c600b0507ea7a75"
  }  
}
