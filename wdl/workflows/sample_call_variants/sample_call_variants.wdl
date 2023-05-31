version 1.0

import "../../tasks/sniffles.wdl" as Sniffles
import "../../tasks/pbsv.wdl" as Pbsv
import "../deepvariant/deepvariant.wdl" as DeepVariant
import "../../tasks/trgt.wdl" as Trgt

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
  if (defined(reference.trgt_tandem_repeat_bed)) {
    call Trgt.trgt {
      input:
        sex = sex,
        bam = aligned_bam.data,
        bam_index = aligned_bam.index,
        reference_fasta = reference.fasta.data,
        reference_index = reference.fasta.index,
        tandem_repeat_bed = select_first([reference.trgt_tandem_repeat_bed]),
        runtime_attributes = default_runtime_attributes
    }
  }

  output {
    File sniffles_snf = sniffles_discover.snf
    Array[File] pbsv_svsigs = pbsv_discover.svsig
    IndexData deepvariant_gvcf = deepvariant.gvcf
    File? trgt_vcf = trgt.repeat_vcf

    # Variant metrics
    # numbers of passing variants by types (hom/het)
    # het/hom ratio by type? 
    # transitions/transversions 
  }
}
