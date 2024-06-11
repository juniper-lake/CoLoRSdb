version 1.0

# Joint call or combine variants for a cohort of samples

import "../tasks/glnexus.wdl" as Glnexus
import "../tasks/sniffles.wdl" as Sniffles
import "../tasks/peddy.wdl" as Peddy
import "../tasks/vcfparser.wdl" as Vcfparser
import "../tasks/bcftools.wdl" as Bcftools
import "../tasks/jasminesv.wdl" as Jasminesv

workflow merge_samples {
  input {
    String cohort_id
    Boolean anonymize_output

    Array[String] sample_ids
    Array[String] sexes
    Array[File] pbsv_vcfs
    Array[File] pbsv_vcf_indexes
    Array[File] unzipped_pbsv_vcfs
    Array[File] deepvariant_gvcfs
    Array[File] deepvariant_gvcf_indexes
    Array[File] sniffles_snfs
    Array[Array[File?]] trgt_vcfs

    ReferenceData reference

    RuntimeAttributes default_runtime_attributes
  }

  scatter (idx in range(length(sample_ids))) {
    String sample_plus_sex = "~{sample_ids[idx]}+~{sexes[idx]}"
  }

  # jasminesv merge pbsv vcfs
  call Jasminesv.jasminesv_merge_svs as jasminesv_merge_pbsv_vcfs {
    input:
      vcfs = unzipped_pbsv_vcfs,
      output_prefix = "~{cohort_id}.~{reference.name}.pbsv.jasminesv",
      runtime_attributes = default_runtime_attributes
  }

  # zip index pbsv/jasminesv
  call Bcftools.reheader_zip_index_jasminesv {
    input:
      vcf = jasminesv_merge_pbsv_vcfs.merged_vcf,
      sample_ids = sample_ids,
      runtime_attributes = default_runtime_attributes
  }

  # pbsv/jasminesv stats per sample
  call Bcftools.structural_variant_stats as pbsv_jasminesv_stats {
    input:
      vcf = reheader_zip_index_jasminesv.zipped_vcf,
      sample_ids = sample_ids,
      autosomes = reference.autosomes,
      runtime_attributes = default_runtime_attributes
  }

  # postprocess pbsv/jasminesv
  call Vcfparser.postprocess_joint_vcf as postprocess_pbsv_vcf {
    input:
      vcf = reheader_zip_index_jasminesv.zipped_vcf,
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sample_plus_sexes = sample_plus_sex,
      non_diploid_regions = reference.non_diploid_regions,
      runtime_attributes = default_runtime_attributes
  }

  # joint call structural variants with sniffles
  call Sniffles.sniffles_call {
    input:
      sample_id = cohort_id,
      snfs = sniffles_snfs,
      reference_name = reference.name,
      reference_fasta = reference.fasta,
      reference_index = reference.fasta_index,
      tr_bed = reference.tandem_repeat_bed,
      runtime_attributes = default_runtime_attributes
  }

  # filter, zip, and index sniffles vcf
  call Bcftools.filter_zip_index_sniffles {
    input:
      vcf = sniffles_call.vcf,
      chromosomes = reference.chromosomes,
      runtime_attributes = default_runtime_attributes
  }

  # sniffles stats per sample
  call Bcftools.structural_variant_stats as sniffles_stats {
    input:
      vcf = filter_zip_index_sniffles.zipped_vcf,
      sample_ids = sample_ids,
      autosomes = reference.autosomes,
      runtime_attributes = default_runtime_attributes
  }

  # postprocess filter/zip/indexed sniffles
  call Vcfparser.postprocess_joint_vcf as postprocess_sniffles_vcf {
    input:
      vcf = filter_zip_index_sniffles.zipped_vcf,
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sample_plus_sexes = sample_plus_sex,
      non_diploid_regions = reference.non_diploid_regions,
      runtime_attributes = default_runtime_attributes
  }

  # joint call small variants with GLnexus
  call Glnexus.glnexus {
    input:
      cohort_id = cohort_id,
      gvcfs = deepvariant_gvcfs,
      reference_name = reference.name,
      runtime_attributes = default_runtime_attributes
  }

  # filter, normalize, and get stats from deepvariant/glnexus vcf
  call Bcftools.filter_norm_deepvariant {
    input:
      vcf = glnexus.vcf,
      vcf_index = glnexus.vcf_index,
      reference_fasta = reference.fasta,
      chromosomes = reference.chromosomes,
      runtime_attributes = default_runtime_attributes
  }

  # postprocess filtered/normalized deepvariant/glnexus vcf
  call Vcfparser.postprocess_joint_vcf as postprocess_deepvariant_vcf {
    input:
      vcf = filter_norm_deepvariant.normalized_vcf,
      cohort_id = cohort_id,
      anonymize_output = anonymize_output,
      sample_plus_sexes = sample_plus_sex,
      non_diploid_regions = reference.non_diploid_regions,
      runtime_attributes = default_runtime_attributes
  }

  # ancestry with peddy from deepvariant vcf
  if (defined(reference.peddy_sites) && defined(reference.peddy_bin)) {
    call Peddy.peddy {
      input:
        cohort_id = cohort_id,
        sample_ids = sample_ids,
        reference_name = reference.name,
        vcf = glnexus.vcf,
        vcf_index = glnexus.vcf_index,
        peddy_sites = select_first([reference.peddy_sites]),
        peddy_bin = select_first([reference.peddy_bin]),
        runtime_attributes = default_runtime_attributes
    }
  }

  # merge trgt vcfs
  if (defined(reference.trgt_tandem_repeat_beds)) {
    Array[Array[File?]] trgt_vcfs_transposed = transpose(trgt_vcfs)

    scatter (idx in range(length(select_first([reference.trgt_tandem_repeat_beds])))) {
      call Vcfparser.merge_trgt_vcfs {
        input:
          trgt_vcfs = select_all(trgt_vcfs_transposed[idx]),
          trgt_bed = select_first([reference.trgt_tandem_repeat_beds])[idx],
          cohort_id = cohort_id,
          reference_name = reference.name,
          reference_index = reference.fasta_index,
          anonymize_output = false,
          runtime_attributes = default_runtime_attributes
      }

      # postprocess trgt if anonymize_output=true, otherwise trgt already fixes ploidy
      if (anonymize_output) {
        call Vcfparser.postprocess_joint_vcf as postprocess_trgt_vcf {
          input:
            vcf = merge_trgt_vcfs.merged_vcf,
            cohort_id = cohort_id,
            anonymize_output = anonymize_output,
            runtime_attributes = default_runtime_attributes
        }
      }
    }
    Array[File] trgt_postprocessed_vcfs_temp = select_all(postprocess_trgt_vcf.postprocessed_vcf)
    Array[File] trgt_postprocessed_vcf_indexes_temp = select_all(postprocess_trgt_vcf.postprocessed_vcf_index)
  }

  output {
    # postprocessed VCFs
    File deepvariant_glnexus_postprocessed_vcf = postprocess_deepvariant_vcf.postprocessed_vcf
    File deepvariant_glnexus_postprocessed_vcf_index = postprocess_deepvariant_vcf.postprocessed_vcf_index
    File pbsv_jasminesv_postprocessed_vcf = postprocess_pbsv_vcf.postprocessed_vcf
    File pbsv_jasminesv_postprocessed_vcf_index = postprocess_pbsv_vcf.postprocessed_vcf_index
    File sniffles_postprocessed_vcf = postprocess_sniffles_vcf.postprocessed_vcf
    File sniffles_postprocessed_vcf_index = postprocess_sniffles_vcf.postprocessed_vcf_index
    Array[File]? trgt_postprocessed_vcfs = trgt_postprocessed_vcfs_temp
    Array[File]? trgt_postprocessed_vcf_indexes = trgt_postprocessed_vcf_indexes_temp

    # original VCFs, so if anonymize_output=false we can access VCFs without ploidy changes
    File deepvariant_glnexus_vcf = filter_norm_deepvariant.normalized_vcf
    File deepvariant_glnexus_vcf_index = filter_norm_deepvariant.normalized_vcf_index
    File sniffles_vcf = filter_zip_index_sniffles.zipped_vcf
    File sniffles_vcf_index = filter_zip_index_sniffles.zipped_vcf_index
    File pbsv_jasminesv_vcf = reheader_zip_index_jasminesv.zipped_vcf
    File pbsv_jasminesv_vcf_index = reheader_zip_index_jasminesv.zipped_vcf_index
    Array[File]? trgt_vcf = merge_trgt_vcfs.merged_vcf
    Array[File]? trgt_vcf_index = merge_trgt_vcfs.merged_vcf_index

    # vcf stats
    File pbsv_jasminesv_vcf_stats = pbsv_jasminesv_stats.stats
    File sniffles_vcf_stats = sniffles_stats.stats

    # ancestry when peddy is run
    File? peddy_het_check = peddy.het_check
    File? peddy_sex_check = peddy.sex_check
    File? peddy_ped_check = peddy.ped_check
    File? peddy_background_pca = peddy.background_pca
    File? peddy_html = peddy.html
    File? peddy_ped = peddy.ped
    File? peddy_vs_html = peddy.vs_html
  }
}
