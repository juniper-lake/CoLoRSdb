version 1.0

# Call variants for a single sample

import "../tasks/sniffles.wdl" as Sniffles
import "../tasks/pbsv.wdl" as Pbsv
import "deepvariant.wdl" as DeepVariant
import "../tasks/trgt.wdl" as Trgt
import "../tasks/bcftools.wdl" as Bcftools

workflow call_variants_by_sample {
  input {
    String sample_id
    String sex
    File aligned_bam
    File aligned_bam_index

    ReferenceData reference

    RuntimeAttributes default_runtime_attributes
  }

  scatter (idx in range(length(reference.chromosomes))) {
    # discover SV signatures with pbsv
    call Pbsv.pbsv_discover {
      input:
        bam = aligned_bam,
        bam_index = aligned_bam_index,
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
        reference_fasta = reference.fasta,
        reference_index = reference.fasta_index,
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
      bam = aligned_bam,
      bam_index = aligned_bam_index,
      reference_fasta = reference.fasta,
      reference_index = reference.fasta_index,
      tandem_repeat_bed = reference.tandem_repeat_bed,
      runtime_attributes = default_runtime_attributes
  }

  # call small variants with deepvariant
  call DeepVariant.deepvariant {
    input:
      sample_id = sample_id,
      aligned_bams = [aligned_bam],
      aligned_bam_indexes = [aligned_bam_index],
      reference_name = reference.name,
      reference_fasta = reference.fasta,
      reference_fasta_index = reference.fasta_index,
      default_runtime_attributes = default_runtime_attributes
  }

  # genotype tandem repeats with trgt
  if (defined(reference.trgt_tandem_repeat_beds)) {
    scatter (trgt_bed in select_first([reference.trgt_tandem_repeat_beds])) {
      call Trgt.trgt {
        input:
          sex = sex,
          bam = aligned_bam,
          bam_index = aligned_bam_index,
          reference_fasta = reference.fasta,
          reference_index = reference.fasta_index,
          tandem_repeat_bed = trgt_bed,
          runtime_attributes = default_runtime_attributes
      }
    }
  }

  output {
    File sniffles_snf = sniffles_discover.snf
    File unzipped_pbsv_vcf = concat_pbsv_vcfs.concatenated_vcf
    File pbsv_vcf = concat_pbsv_vcfs.concatenated_zipped_vcf
    File pbsv_vcf_index = concat_pbsv_vcfs.concatenated_zipped_vcf_index
    File deepvariant_gvcf = deepvariant.gvcf
    File deepvariant_gvcf_index = deepvariant.gvcf_index
    Array[File]? trgt_vcf = trgt.repeat_vcf
  }
}
