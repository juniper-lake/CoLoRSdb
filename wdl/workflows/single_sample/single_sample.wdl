version 1.0

import "../../tasks/pbmm2.wdl" as pbmm2
import "../../tasks/somalier.wdl" as somalier
import "../../tasks/sniffles.wdl" as sniffles
import "../../tasks/pbsv.wdl" as pbsv
import "../deepvariant/deepvariant.wdl" as DeepVariant
import "../../tasks/utils.wdl" as utils
import "../../tasks/trgt.wdl" as Trgt

workflow sample_analysis {

  input {
    Sample sample
    ReferenceData reference
    Float min_relatedness_sample_swap
    String deepvariant_version
    RuntimeAttributes default_runtime_attributes
  }

   # align each movie with pbmm2
  scatter (idx in range(length(sample.movies))) {     
    call pbmm2.pbmm2_align {
        input: 
          movie = sample.movies[idx],
          sample_id = sample.sample_id,
          reference_name = reference.name,
          reference_fasta = reference.fasta.data,
          reference_index = reference.fasta.index,
          runtime_attributes = default_runtime_attributes
    }
    
    IndexData aligned_bam = {
			"data": pbmm2_align.bam,
			"index": pbmm2_align.bam_index
		}
  }

  # check for sample swaps and confirm gender of sample
  call somalier.somalier_relate {
    input:
      sample_id = sample.sample_id,
      bams = pbmm2_align.bam,
      bam_indexes = pbmm2_align.bam_index,
      reference_fasta = reference.fasta.data,
      reference_index = reference.fasta.index,
      sample_swap_sites_vcf = reference.somalier_sites_vcf,
      runtime_attributes = default_runtime_attributes
  }

  if (somalier_relate.min_relatedness >= min_relatedness_sample_swap) {
    
    if (somalier_relate.inferred_sex != "U") {

      scatter (idx in range(length(pbmm2_align.bam))) {
        # call structural variants with sniffles
        call sniffles.sniffles_discover {
          input:
            sample_id = sample.sample_id,
            bam = pbmm2_align.bam[idx],
            bam_index = pbmm2_align.bam_index[idx],
            reference_fasta = reference.fasta.data,
            reference_index = reference.fasta.index,
            tr_bed = reference.tandem_repeat_bed,
            runtime_attributes = default_runtime_attributes
        }
      }

      scatter (idx in range(length(pbmm2_align.bam))) {
        # discover SV signatures with pbsv
        call pbsv.pbsv_discover {
          input:
            bam = pbmm2_align.bam[idx],
            bam_index = pbmm2_align.bam_index[idx],
            tr_bed = reference.tandem_repeat_bed,
            runtime_attributes = default_runtime_attributes
        }
      }
      Array[File] svsig = flatten(pbsv_discover.svsigs)

      call DeepVariant.deepvariant {
        input:
          sample_id = sample.sample_id,
          aligned_bams = aligned_bam,
          reference_name = reference.name,
          reference_fasta = reference.fasta,
          deepvariant_version = deepvariant_version,
          default_runtime_attributes = default_runtime_attributes
      }

      call utils.merge_bams {
        input:
          bams = pbmm2_align.bam,
          output_bam_name = "~{sample.sample_id}.~{reference.name}.bam",
          runtime_attributes = default_runtime_attributes
      }

      IndexData merged_bam = {
        "data": merge_bams.merged_bam, 
        "index": merge_bams.merged_bam_index
      }

      call Trgt.trgt {
        input:
          bam = merge_bams.merged_bam,
          bam_index = merge_bams.merged_bam_index,
          reference_fasta = reference.fasta.data,
          reference_index = reference.fasta.index,
          tandem_repeat_bed = reference.tandem_repeat_bed,
          runtime_attributes = default_runtime_attributes
      }

      IndexData trgt_repeat_vcf = {
        "data": trgt.repeat_vcf,
        "index": trgt.repeat_vcf_index
      }

      # keep track of sample IDs that pass QC
      String qc_pass = sample.sample_id
    }
  }

  if (somalier_relate.min_relatedness < min_relatedness_sample_swap) {
    # keep track of sample IDs that are excluded for sample swap between movies
    String qc_fail_relatedness = sample.sample_id
  }
  if (somalier_relate.inferred_sex == "U") {
    # keep track of sample IDs that are excluded for unknown sex
    String qc_fail_sex = sample.sample_id
  }

  output {
    # QC pass
    IndexData? pbmm2_merged_bam = merged_bam
    Array[File]? sniffles_snfs = sniffles_discover.snf
    Array[File]? pbsv_svsigs = svsig
    IndexData? deepvariant_gvcf = deepvariant.gvcf
    IndexData? trgt_vcf = trgt_repeat_vcf
    String? qc_pass_id = qc_pass
    
    # QC fail
    String? qc_fail_relatedness_id = qc_fail_relatedness
    String? qc_fail_sex_id = qc_fail_sex
  }
}
