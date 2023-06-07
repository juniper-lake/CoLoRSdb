version 1.0

import "../../tasks/glnexus.wdl" as Glnexus
import "../../tasks/pbsv.wdl" as Pbsv
import "../../tasks/sniffles.wdl" as Sniffles
import "../../tasks/hificnv.wdl" as Hificnv
import "../../tasks/hiphase.wdl" as Hiphase
import "../../tasks/peddy.wdl" as Peddy
import "../../tasks/vcfparser.wdl" as Vcfparser
import "../../tasks/utils.wdl" as Utils

workflow cohort_combine_samples {
  input {
    String cohort_id
    Boolean anonymize_output

    Array[String] sample_ids
    Array[String] sexes
    Array[IndexData] aligned_bams
    Array[Array[File]] svsigs
    Array[IndexData] gvcfs
    Array[File] snfs
    Array[File?] trgt_vcfs
    Array[File?] hificnv_vcfs

    ReferenceData reference

    RuntimeAttributes default_runtime_attributes 
  }

  scatter (gvcf_object in gvcfs) {
		File gvcf = gvcf_object.data
		File gvcf_index = gvcf_object.index
	}

  # joint call small variants with GLnexus
  call Glnexus.glnexus {
    input:
      cohort_id = cohort_id,
      gvcfs = gvcf,
      gvcf_indexes = gvcf_index,
      reference_name = reference.name,
      runtime_attributes = default_runtime_attributes
  }

  # ancestry with peddy
  if (defined(reference.peddy_sites) && defined(reference.peddy_bin)) {
    call Peddy.peddy {
      input:
        cohort_id = cohort_id,
        sample_ids = sample_ids,
        vcf = glnexus.vcf,
        vcf_index = glnexus.vcf_index,
        peddy_sites = select_first([reference.peddy_sites]),
        peddy_bin = select_first([reference.peddy_bin]),
        runtime_attributes = default_runtime_attributes
    }
  }

  # joint call structural variants with pbsv
  scatter (idx in range(length(reference.chromosomes))) {
    call Pbsv.pbsv_call {
      input:
        sample_id = cohort_id,
        svsigs = svsigs[idx],
        sample_count = length(sample_ids),
        region = reference.chromosomes[idx],
        reference_name = reference.name,
        reference_fasta = reference.fasta.data,
        reference_index = reference.fasta.index,
        runtime_attributes = default_runtime_attributes
    }
  }

  call Utils.concat_vcfs {
    input:
      vcfs = pbsv_call.vcf,
      output_vcf_name = "~{cohort_id}.~{reference.name}.pbsv.vcf",
      runtime_attributes = default_runtime_attributes
  }

  # joint call structural variants with sniffles
  call Sniffles.sniffles_call {
    input:
      sample_id = cohort_id,
      snfs = snfs,
      reference_name = reference.name,
      reference_fasta = reference.fasta.data,
      reference_index = reference.fasta.index,
      tr_bed = reference.tandem_repeat_bed,
      runtime_attributes = default_runtime_attributes
  }

  scatter (bam_object in aligned_bams) {
		File aligned_bam = bam_object.data
		File aligned_bam_index = bam_object.index
	}

  if (!anonymize_output) {
    # phase pbsv and deepvariant vcfs
    call Hiphase.hiphase {
      input:
        cohort_id = cohort_id,
        sample_ids = sample_ids,
        aligned_bams = aligned_bam,
        aligned_bam_indexes = aligned_bam_index,
        deepvariant_vcf = glnexus.vcf,
        deepvariant_vcf_index = glnexus.vcf_index,
        pbsv_vcf = concat_vcfs.concatenated_vcf,
        reference_fasta = reference.fasta.data,
        reference_index = reference.fasta.index,
        runtime_attributes = default_runtime_attributes
    }
  }

  scatter (idx in range(length(sample_ids))) {
    String sample_plus_sex = "~{sample_ids[idx]}+~{sexes[idx]}}"
  }

  # pbsv
  call Vcfparser.postprocess_joint_vcf as postprocess_pbsv_vcf {
    input:
      vcf = select_first([hiphase.pbsv_output_vcf, concat_vcfs.concatenated_vcf]),
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sample_plus_sexes = sample_plus_sex,
      haploid_bed = reference.haploid_bed,
      runtime_attributes = default_runtime_attributes
  }

  # deepvariant
  call Vcfparser.postprocess_joint_vcf as postprocess_deepvariant_vcf {
    input:
      vcf = select_first([hiphase.deepvariant_output_vcf, glnexus.vcf]),
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sample_plus_sexes = sample_plus_sex,
      haploid_bed = reference.haploid_bed,
      runtime_attributes = default_runtime_attributes
  }

  # sniffles
  call Vcfparser.postprocess_joint_vcf as postprocess_sniffles_vcf {
    input:
      vcf = sniffles_call.vcf,
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sample_plus_sexes = sample_plus_sex,
      haploid_bed = reference.haploid_bed,
      runtime_attributes = default_runtime_attributes
  }

  # trgt 
  if (length(trgt_vcfs) > 0) {
    call Vcfparser.merge_trgt_vcfs {
      input:
        trgt_vcfs = select_all(trgt_vcfs),
        cohort_id = cohort_id,
        reference_name = reference.name,
        anonymize_output = anonymize_output,
        runtime_attributes = default_runtime_attributes
    }

    IndexData cohort_trgt = { 
      "data": merge_trgt_vcfs.merged_vcf,
      "index": merge_trgt_vcfs.merged_vcf_index
      }
  }
  
  # hificnv
  if (length(hificnv_vcfs) > 0) {
    call Hificnv.merge_hificnv_vcfs {
      input:
        cnv_vcfs = select_all(hificnv_vcfs),
        cohort_id = cohort_id,
        reference_name = reference.name,
        runtime_attributes = default_runtime_attributes
    }

    call Vcfparser.postprocess_joint_vcf as postprocess_hificnv_vcf {
      input:
        vcf = merge_hificnv_vcfs.merged_cnv_vcf,
        cohort_id = cohort_id,
        anonymize_output = anonymize_output,
        runtime_attributes = default_runtime_attributes
    }

    IndexData cohort_hificnv = {
      "data": postprocess_hificnv_vcf.postprocessed_vcf,
      "index": postprocess_hificnv_vcf.postprocessed_vcf_index
    }
  }


  output {
    IndexData cohort_deepvariant_vcf = { 
      "data": postprocess_deepvariant_vcf.postprocessed_vcf,
      "index": postprocess_deepvariant_vcf.postprocessed_vcf_index
      }
    IndexData cohort_pbsv_vcf = { 
      "data": postprocess_pbsv_vcf.postprocessed_vcf,
      "index": postprocess_pbsv_vcf.postprocessed_vcf_index 
      }
    IndexData cohort_sniffles_vcf = { 
      "data": postprocess_sniffles_vcf.postprocessed_vcf,
      "index": postprocess_sniffles_vcf.postprocessed_vcf_index
      }
    IndexData? cohort_trgt_vcf = cohort_trgt
    IndexData? cohort_hificnv_vcf = cohort_hificnv
    Array[File] pbsv_call_logs = pbsv_call.log # for testing memory usage
    # File variant_summary

    File? peddy_het_check = peddy.het_check
    File? peddy_sex_check = peddy.sex_check
    File? peddy_ped_check = peddy.ped_check
    File? peddy_background_pca = peddy.background_pca
    File? peddy_html = peddy.html
    File? peddy_ped = peddy.ped
    File? peddy_vs_html = peddy.vs_html

    # Aggregate = false
    File? hiphase_stats = hiphase.stats
    File? hiphase_blocks = hiphase.blocks
    File? hiphase_summary = hiphase.summary
  }
}
