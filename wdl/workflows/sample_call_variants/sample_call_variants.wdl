version 1.0

import "../../tasks/pbmm2.wdl" as Pbmm2
import "../../tasks/somalier.wdl" as Somalier
import "../../tasks/sniffles.wdl" as Sniffles
import "../../tasks/pbsv.wdl" as Pbsv
import "../sample_deepvariant/sample_deepvariant.wdl" as DeepVariant
import "../../tasks/utils.wdl" as Utils
import "../../tasks/trgt.wdl" as Trgt

workflow sample_call_variants {

  input {
    String sample_id
    IndexData aligned_bam

    ReferenceData reference

    String deepvariant_version

    RuntimeAttributes default_runtime_attributes
  }

  scatter (idx in range(length(reference.chromosomes))) {
    # discover SV signatures with pbsv
    call Pbsv.pbsv_discover {
      input:
        bam = merged_bam.data,
        bam_index = merged_bam.index,
        region = reference.chromosomes[idx],
        tandem_repeat_bed = reference.tandem_repeat_bed,
        runtime_attributes = default_runtime_attributes
    }
  }

  # call structural variants with sniffles
  call Sniffles.sniffles_discover {
    input:
      sample_id = sample.sample_id,
      bam = merged_bam.data,
      bam_index = merged_bam.index,
      reference_fasta = reference.fasta.data,
      reference_index = reference.fasta.index,
      tr_bed = reference.tandem_repeat_bed,
      runtime_attributes = default_runtime_attributes
  }

  # call small variants with deepvariant
  call DeepVariant.deepvariant {
    input:
      sample_id = sample.sample_id,
      aligned_bams = aligned_bam,
      reference_name = reference.name,
      reference_fasta = reference.fasta,
      deepvariant_version = deepvariant_version,
      default_runtime_attributes = default_runtime_attributes
  }

  # genotype tandem repeats with trgt
  call Trgt.trgt {
    input:
      sex = somalier_relate.inferred_sex,
      bam = merge_bams.merged_bam,
      bam_index = merge_bams.merged_bam_index,
      reference_fasta = reference.fasta.data,
      reference_index = reference.fasta.index,
      tandem_repeat_bed = reference.trgt_tandem_repeat_bed,
      runtime_attributes = default_runtime_attributes
  }

  output {
    Array[File] sniffles_snfs = sniffles_discover.snf
    Array[File] pbsv_svsigs = pbsv_discover.svsig
    IndexData deepvariant_gvcf = deepvariant.gvcf
    IndexData trgt_vcf = { "data": trgt.repeat_vcf, "index": trgt.repeat_vcf_index }

    # Variant metrics
    # numbers of passing variants by types (hom/het)
    # het/hom ratio by type? 
    # transitions/transversions 
  }
}
