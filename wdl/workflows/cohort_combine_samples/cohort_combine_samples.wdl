version 1.0

# Joint call or combine variants for a cohort of samples

import "../../tasks/glnexus.wdl" as Glnexus
import "../../tasks/sniffles.wdl" as Sniffles
import "../../tasks/hificnv.wdl" as Hificnv
import "../../tasks/peddy.wdl" as Peddy
import "../../tasks/vcfparser.wdl" as Vcfparser
import "../../tasks/bcftools.wdl" as Bcftools
import "../../tasks/jasminesv.wdl" as Jasminesv

workflow cohort_combine_samples {
  input {
    String cohort_id
    Boolean anonymize_output

    Array[String] sample_ids
    Array[String] sexes
    Array[IndexData] pbsv_vcfs
    Array[File] unzipped_pbsv_vcfs
    Array[IndexData] deepvariant_gvcfs
    Array[File] sniffles_snfs
    Array[Array[File?]] trgt_vcfs
    Array[IndexData?] hificnv_vcfs

    ReferenceData reference

    RuntimeAttributes default_runtime_attributes
  }

  scatter (gvcf_object in deepvariant_gvcfs) {
    File deepvariant_gvcf = gvcf_object.data
    File deepvariant_gvcf_index = gvcf_object.index
  }

  scatter (vcf_object in pbsv_vcfs) {
    File pbsv_vcf = vcf_object.data
    File pbsv_vcf_index = vcf_object.index
  }

  scatter (idx in range(length(sample_ids))) {
    String sample_plus_sex = "~{sample_ids[idx]}+~{sexes[idx]}"
  }

  # strict merge pbsv vcfs with bcftools
  call Bcftools.merge_vcfs as strict_merge_pbsv_vcfs {
    input:
      vcfs = pbsv_vcf,
      vcf_indexes = pbsv_vcf_index,
      output_prefix = "~{cohort_id}.~{reference.name}.pbsv.strictmerge",
      runtime_attributes = default_runtime_attributes
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
      reference_fasta = reference.fasta.data,
      reference_index = reference.fasta.index,
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
      gvcfs = deepvariant_gvcf,
      gvcf_indexes = deepvariant_gvcf_index,
      reference_name = reference.name,
      runtime_attributes = default_runtime_attributes
  }

  # filter, normalize, and get stats from deepvariant/glnexus vcf
  call Bcftools.filter_norm_deepvariant {
    input:
      vcf = glnexus.vcf,
      vcf_index = glnexus.vcf_index,
      reference_fasta = reference.fasta.data,
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
          trgt_bed = reference.trgt_tandem_repeat_beds[idx],
          cohort_id = cohort_id,
          reference_name = reference.name,
          reference_index = reference.fasta.index,
          anonymize_output = false,
          runtime_attributes = default_runtime_attributes
      }

      IndexData merged_trgt = {
        "data": merge_trgt_vcfs.merged_vcf,
        "index": merge_trgt_vcfs.merged_vcf_index
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

        IndexData postprocessed_trgt = {
          "data": postprocess_trgt_vcf.postprocessed_vcf,
          "index": postprocess_trgt_vcf.postprocessed_vcf_index
        }
      }
    }

    Array[IndexData] postprocessed_trgts = select_all(postprocessed_trgt)
  }

  if (length(hificnv_vcfs) > 0) {
      scatter (vcf_object in select_all(hificnv_vcfs)) {
        File hificnv_vcf_data = vcf_object.data
        File hificnv_vcf_index = vcf_object.index
      }

    # merge hificnv vcfs
    call Hificnv.merge_hificnv_vcfs {
      input:
        cnv_vcfs = hificnv_vcf_data,
        cnv_vcf_indexes = hificnv_vcf_index,
        cohort_id = cohort_id,
        reference_name = reference.name,
        runtime_attributes = default_runtime_attributes
    }

    IndexData merged_hificnv = {
      "data": merge_hificnv_vcfs.merged_cnv_vcf,
      "index": merge_hificnv_vcfs.merged_cnv_vcf_index
    }

    # postprocess hificnv
    call Vcfparser.postprocess_joint_vcf as postprocess_hificnv_vcf {
      input:
        vcf = merge_hificnv_vcfs.merged_cnv_vcf,
        cohort_id = cohort_id,
        sample_plus_sexes = sample_plus_sex,
        non_diploid_regions = reference.non_diploid_regions,
        anonymize_output = anonymize_output,
        runtime_attributes = default_runtime_attributes
    }

    IndexData postprocessed_hificnv = {
      "data": postprocess_hificnv_vcf.postprocessed_vcf,
      "index": postprocess_hificnv_vcf.postprocessed_vcf_index
    }
  }

  output {
    # postprocessed VCFs
    IndexData deepvariant_glnexus_postprocessed_vcf = {
      "data": postprocess_deepvariant_vcf.postprocessed_vcf,
      "index": postprocess_deepvariant_vcf.postprocessed_vcf_index
    }
    IndexData pbsv_jasminesv_postprocessed_vcf = {
      "data": postprocess_pbsv_vcf.postprocessed_vcf,
      "index": postprocess_pbsv_vcf.postprocessed_vcf_index
    }
    IndexData sniffles_postprocessed_vcf = {
      "data": postprocess_sniffles_vcf.postprocessed_vcf,
      "index": postprocess_sniffles_vcf.postprocessed_vcf_index
    }
    Array[IndexData]? trgt_postprocessed_vcf = postprocessed_trgts
    IndexData? hificnv_postprocessed_vcf = postprocessed_hificnv

    # original VCFs, so if anonymize_output=false we can access VCFs without ploidy changes
    IndexData deepvariant_glnexus_vcf = {
      "data": filter_norm_deepvariant.normalized_vcf,
      "index": filter_norm_deepvariant.normalized_vcf_index
    }
    IndexData sniffles_vcf = {
      "data": filter_zip_index_sniffles.zipped_vcf,
      "index": filter_zip_index_sniffles.zipped_vcf_index
    }
    IndexData pbsv_strictmerge_vcf = {
      "data": strict_merge_pbsv_vcfs.merged_vcf,
      "index": strict_merge_pbsv_vcfs.merged_vcf_index
    }
    IndexData pbsv_jasminesv_vcf = {
      "data": reheader_zip_index_jasminesv.zipped_vcf,
      "index": reheader_zip_index_jasminesv.zipped_vcf_index
    }
    Array[IndexData]? trgt_vcf = merged_trgt
    IndexData? hificnv_vcf = merged_hificnv

    # VCF stats
    File pbsv_jasminesv_vcf_stats = pbsv_jasminesv_stats.stats
    File sniffles_vcf_stats = sniffles_stats.stats

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
