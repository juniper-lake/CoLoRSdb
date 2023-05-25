version 1.0

import "../../tasks/glnexus.wdl" as Glnexus
import "../../tasks/pbsv.wdl" as Pbsv
import "../../tasks/sniffles.wdl" as Sniffles
import "../../tasks/hiphase.wdl" as Hiphase
import "../../tasks/vcfparser.wdl" as VcfParser
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

  # pbsv
  call VcfParser.postprocess_joint_vcf as postprocess_pbsv_vcf {
    input:
      vcf = select_first([hiphase.pbsv_output_vcf, concat_vcfs.concatenated_vcf]),
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sexes = sexes,
      runtime_attributes = default_runtime_attributes
  }

  # deepvariant
  call VcfParser.postprocess_joint_vcf as postprocess_deepvariant_vcf {
    input:
      vcf = select_first([hiphase.deepvariant_output_vcf, glnexus.vcf]),
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sexes = sexes,
      runtime_attributes = default_runtime_attributes
  }

  # sniffles
  call VcfParser.postprocess_joint_vcf as postprocess_sniffles_vcf {
    input:
      vcf = sniffles_call.vcf,
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sexes = sexes,
      runtime_attributes = default_runtime_attributes
  }

  # trgt 
  if (length(trgt_vcfs) > 0) {
    call VcfParser.merge_trgt_vcfs {
      input:
        trgt_vcfs = select_all(trgt_vcfs),
        cohort_id = cohort_id,
        reference_name = reference.name,
        anonymize_output = anonymize_output,
        runtime_attributes = default_runtime_attributes
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
    IndexData? cohort_trgt_vcf = { 
      "data": select_first([merge_trgt_vcfs.merged_trgt_vcf]),
      "index": select_first([merge_trgt_vcfs.merge_trgt_vcf_index])
      }
    Array[File] pbsv_call_logs = pbsv_call.log # for testing memory usage
    # File variant_summary

    # Aggregate = false
    File? hiphase_stats = hiphase.stats
    File? hiphase_blocks = hiphase.blocks
    File? hiphase_summary = hiphase.summary
  }
    
}



