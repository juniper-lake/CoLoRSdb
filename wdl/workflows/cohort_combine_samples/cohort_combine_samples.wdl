version 1.0

# Joint call or combine variants for a cohort of samples

import "../../tasks/glnexus.wdl" as Glnexus
import "../../tasks/pbsv.wdl" as Pbsv
import "../../tasks/sniffles.wdl" as Sniffles
import "../../tasks/hificnv.wdl" as Hificnv
import "../../tasks/peddy.wdl" as Peddy
import "../../tasks/vcfparser.wdl" as Vcfparser
import "../../tasks/bcftools.wdl" as Bcftools

workflow cohort_combine_samples {
  input {
    String cohort_id
    Boolean anonymize_output
    Int max_samples_pbsv_call

    Array[String] sample_ids
    Array[String] sexes
    Array[Array[File]] pbsv_svsigs
    Array[IndexData] deepvariant_gvcfs
    Array[File] sniffles_snfs
    Array[File?] trgt_vcfs
    Array[IndexData?] hificnv_vcfs

    ReferenceData reference

    Int? pbsv_call_mem_gb
    Int? glnexus_mem_gb
    Int? sniffles_call_mem_gb

    RuntimeAttributes default_runtime_attributes 
  }

  scatter (gvcf_object in deepvariant_gvcfs) {
    File deepvariant_gvcf = gvcf_object.data
    File deepvariant_gvcf_index = gvcf_object.index
  }

  # joint call small variants with GLnexus
  call Glnexus.glnexus {
    input:
      cohort_id = cohort_id,
      gvcfs = deepvariant_gvcf,
      gvcf_indexes = deepvariant_gvcf_index,
      reference_name = reference.name,
      mem_gb = glnexus_mem_gb,
      runtime_attributes = default_runtime_attributes
  }

  # get stats from deepvariant vcf
  call Bcftools.small_variant_stats {
    input:
      vcf = glnexus.vcf,
      sample_ids = sample_ids,
      reference = reference.fasta.data,
      autosomes = reference.autosomes,
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

  scatter (idx in range(length(sample_ids))) {
    String sample_plus_sex = "~{sample_ids[idx]}+~{sexes[idx]}"
  }

  # limit # of samples being routed to pbsv because it will fail with too many samples
  # get max size of samples for pbsv to process
  # the following commented line didn't work (always rounded down despite ceil function)
  # Int n_pbsv_call_groups  = ceil(length(sample_ids)/max_samples_pbsv_call)
  Int n_pbsv_call_groups  = if length(sample_ids) % max_samples_pbsv_call > 0 then floor(length(sample_ids)/max_samples_pbsv_call) + 1 else length(sample_ids)/max_samples_pbsv_call

  scatter (group_idx in range(n_pbsv_call_groups)) {
    # subsample based on sample index
    scatter (sample_idx in range(length(sample_ids))) {
      Int sample_group = sample_idx % n_pbsv_call_groups
      if (sample_group == group_idx) {
        Array[File] group_svsigs = pbsv_svsigs[sample_idx]
        String group_sample_ids = sample_ids[sample_idx]
      }
    }

    # transpose svsigs to order by chromosome instead of sample
    Array[Array[File]] group_svsigs_transposed = transpose(select_all(group_svsigs))

    # joint call structural variants with pbsv
    scatter (idx in range(length(reference.chromosomes))) {
      call Pbsv.pbsv_call {
        input:
          sample_id = "~{cohort_id}.group_~{group_idx}",
          svsigs = group_svsigs_transposed[idx],
          sample_count = length(group_sample_ids),
          region = reference.chromosomes[idx],
          reference_name = reference.name,
          reference_fasta = reference.fasta.data,
          reference_index = reference.fasta.index,
          mem_gb = pbsv_call_mem_gb,
          runtime_attributes = default_runtime_attributes
      }
    }

    # concat chromosome-specific pbsv calls into single vcf
    call Bcftools.concat_vcfs {
      input:
        vcfs = pbsv_call.vcf,
        output_vcf_name = "~{cohort_id}.group_~{group_idx}.~{reference.name}.pbsv.vcf",
        runtime_attributes = default_runtime_attributes
    }

    # pbsv stats per sample
    call Bcftools.structural_variant_stats as pbsv_stats {
      input:
        vcf = concat_vcfs.concatenated_vcf,
        sample_ids = sample_ids,
        autosomes = reference.autosomes,
        runtime_attributes = default_runtime_attributes
    }

    # postprocess pbsv
    call Vcfparser.postprocess_joint_vcf as postprocess_pbsv_vcf {
      input:
        vcf = concat_vcfs.concatenated_vcf,
        cohort_id = cohort_id,
        anonymize_output = anonymize_output,
        sample_plus_sexes = sample_plus_sex,
        non_diploid_regions = reference.non_diploid_regions,
        runtime_attributes = default_runtime_attributes
    }
  
    IndexData pbsv_vcf = { 
      "data": postprocess_pbsv_vcf.postprocessed_vcf,
      "index": postprocess_pbsv_vcf.postprocessed_vcf_index 
    }
  }

  # joint call structural variants with sniffles
  call Sniffles.sniffles_call {
    input:
      sample_id = cohort_id,
      snfs = sniffles_snfs,
      reference_name = reference.name,
      reference_fasta = reference.fasta.data,
      reference_index = reference.fasta.index,
      tr_bed = reference.tandem_repeat_bed,
      mem_gb = sniffles_call_mem_gb,
      runtime_attributes = default_runtime_attributes
  }

  call Bcftools.structural_variant_stats as sniffles_stats {
    input:
      vcf = sniffles_call.vcf,
      sample_ids = sample_ids,
      autosomes = reference.autosomes,
      runtime_attributes = default_runtime_attributes
  }

  # postprocess deepvariant
  call Vcfparser.postprocess_joint_vcf as postprocess_deepvariant_vcf {
    input:
      vcf = glnexus.vcf,
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sample_plus_sexes = sample_plus_sex,
      non_diploid_regions = reference.non_diploid_regions,
      runtime_attributes = default_runtime_attributes
  }

  # sniffles
  call Vcfparser.postprocess_joint_vcf as postprocess_sniffles_vcf {
    input:
      vcf = sniffles_call.vcf,
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sample_plus_sexes = sample_plus_sex,
      non_diploid_regions = reference.non_diploid_regions,
      runtime_attributes = default_runtime_attributes
  }

  # postprocess trgt 
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
  
  # postprocess hificnv
  if (length(hificnv_vcfs) > 0) {
      scatter (vcf_object in select_all(hificnv_vcfs)) {
        File hificnv_vcf = vcf_object.data
        File hificnv_vcf_index = vcf_object.index
      }

    call Hificnv.merge_hificnv_vcfs {
      input:
        cnv_vcfs = hificnv_vcf,
        cnv_vcf_indexes = hificnv_vcf_index,
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
    # VCFs
    IndexData cohort_deepvariant_vcf = { 
      "data": postprocess_deepvariant_vcf.postprocessed_vcf,
      "index": postprocess_deepvariant_vcf.postprocessed_vcf_index
      }
    Array[IndexData]+ cohort_pbsv_vcfs = pbsv_vcf
    IndexData cohort_sniffles_vcf = { 
      "data": postprocess_sniffles_vcf.postprocessed_vcf,
      "index": postprocess_sniffles_vcf.postprocessed_vcf_index
      }
    IndexData? cohort_trgt_vcf = cohort_trgt
    IndexData? cohort_hificnv_vcf = cohort_hificnv

    # Logs
    Array[Array[File]]+ pbsv_call_logs = pbsv_call.log # for testing memory usage
    
    # VCF stats
    File cohort_deepvariant_vcf_stats = small_variant_stats.stats
    Array[File]+ cohort_pbsv_vcf_stats = pbsv_stats.stats
    File cohort_sniffles_vcf_stats = sniffles_stats.stats
    
    # Ancestry when peddy is run
    File? peddy_het_check = peddy.het_check
    File? peddy_sex_check = peddy.sex_check
    File? peddy_ped_check = peddy.ped_check
    File? peddy_background_pca = peddy.background_pca
    File? peddy_html = peddy.html
    File? peddy_ped = peddy.ped
    File? peddy_vs_html = peddy.vs_html
  }
}
