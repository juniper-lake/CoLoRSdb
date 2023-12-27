version 1.0

# Call variants for a single sample

import "../../tasks/sniffles.wdl" as Sniffles
import "../../tasks/pbsv.wdl" as Pbsv
import "../deepvariant/deepvariant.wdl" as DeepVariant
import "../../tasks/trgt.wdl" as Trgt
import "../../tasks/hificnv.wdl" as Hificnv
import "../../tasks/bcftools.wdl" as Bcftools

workflow sample_call_variants {
  input {
    String sample_id
    IndexData aligned_bam
    String sex

    ReferenceData reference

    String deepvariant_version

    RuntimeAttributes default_runtime_attributes
  }

  scatter (idx in range(length(reference.chromosomes))) {
    # discover SV signatures with pbsv
    call Pbsv.pbsv_discover {
      input:
        bam = aligned_bam.data,
        bam_index = aligned_bam.index,
        region = reference.chromosomes[idx],
        tandem_repeat_bed = reference.tandem_repeat_bed,
        runtime_attributes = default_runtime_attributes
    }

    # call SVs from svsigs with pbsv
    call Pbsv.pbsv_call {
      input:
        sample_id = sample_id,
        svsigs = [pbsv_discover.svsig],
        region = reference.chromosomes[idx],
        reference_name = reference.name,
        reference_fasta = reference.fasta.data,
        reference_index = reference.fasta.index,
        runtime_attributes = default_runtime_attributes
    }
  }

  # concat chromosome-specific pbsv calls into single vcf
  call Bcftools.concat_pbsv_vcfs {
    input:
      vcfs = pbsv_call.vcf,
      output_prefix = "~{sample_id}.~{reference.name}.pbsv",
      runtime_attributes = default_runtime_attributes
  }

  # call structural variants with sniffles
  call Sniffles.sniffles_discover {
    input:
      sample_id = sample_id,
      bam = aligned_bam.data,
      bam_index = aligned_bam.index,
      reference_fasta = reference.fasta.data,
      reference_index = reference.fasta.index,
      tandem_repeat_bed = reference.tandem_repeat_bed,
      runtime_attributes = default_runtime_attributes
  }

  # call small variants with deepvariant
  call DeepVariant.deepvariant {
    input:
      sample_id = sample_id,
      aligned_bams = [aligned_bam],
      reference_name = reference.name,
      reference_fasta = reference.fasta,
      deepvariant_version = deepvariant_version,
      default_runtime_attributes = default_runtime_attributes
  }

  # genotype tandem repeats with trgt
  if (defined(reference.trgt_tandem_repeat_beds)) {
    scatter (trgt_bed in select_first([reference.trgt_tandem_repeat_beds])) {
      call Trgt.trgt {
        input:
          sex = sex,
          bam = aligned_bam.data,
          bam_index = aligned_bam.index,
          reference_fasta = reference.fasta.data,
          reference_index = reference.fasta.index,
          tandem_repeat_bed = trgt_bed,
          runtime_attributes = default_runtime_attributes
      }
    }
  }

  if (defined(reference.hificnv_exclude_bed)
    && defined(reference.hificnv_expected_bed_male)
    && defined(reference.hificnv_expected_bed_female)) {
      call Hificnv.hificnv {
        input:
          sample_id = sample_id,
          sex = sex,
          bam = aligned_bam.data,
          bam_index = aligned_bam.index,
          small_variant_vcf = deepvariant.vcf.data,
          small_variant_vcf_index = deepvariant.vcf.index,
          reference_name = reference.name,
          reference = reference.fasta.data,
          reference_index = reference.fasta.index,
          exclude_bed = select_first([reference.hificnv_exclude_bed]).data,
          exclude_bed_index = select_first([reference.hificnv_exclude_bed]).index,
          expected_bed_male = select_first([reference.hificnv_expected_bed_male]),
          expected_bed_female = select_first([reference.hificnv_expected_bed_female]),
          runtime_attributes = default_runtime_attributes
      }

      IndexData hificnv_indexed_vcf = {
        "data": hificnv.cnv_vcf,
        "index": hificnv.cnv_vcf_index}
    }

  output {
    File sniffles_snf = sniffles_discover.snf
    File unzipped_pbsv_vcf = concat_pbsv_vcfs.concatenated_vcf
    IndexData pbsv_vcf = {
      "data": concat_pbsv_vcfs.concatenated_zipped_vcf,
      "index": concat_pbsv_vcfs.concatenated_zipped_vcf_index
      }
    IndexData deepvariant_gvcf = deepvariant.gvcf
    Array[File]? trgt_vcf = trgt.repeat_vcf
    IndexData? hificnv_vcf = hificnv_indexed_vcf
  }
}
