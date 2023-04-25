version 1.0

import "../../tasks/pbmm2.wdl" as Pbmm2
import "../../tasks/somalier.wdl" as Somalier
import "../../tasks/sniffles.wdl" as Sniffles
import "../../tasks/pbsv.wdl" as Pbsv
import "../deepvariant/deepvariant.wdl" as DeepVariant
import "../../tasks/utils.wdl" as Utils
import "../../tasks/trgt.wdl" as Trgt

workflow single_sample {

  input {
    Sample sample
    ReferenceData reference
    Float min_relatedness_sample_swap
    String deepvariant_version
    RuntimeAttributes default_runtime_attributes
  }

  # align each movie with pbmm2
  scatter (idx in range(length(sample.movies))) {     
    call Pbmm2.pbmm2_align {
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

    # get movie name from movie path
    String movie_name = sub(basename(sample.movies[idx]), "\\..*", "")

    # extract somalier sites from aligned bam
    call Somalier.somalier_extract {
      input:
        sample_id = sample.sample_id,
        movie_name = movie_name,
        bam = aligned_bam.data,
        bam_index = aligned_bam.index,
        reference_fasta = reference.fasta.data,
        reference_index = reference.fasta.index,
        somalier_sites_vcf = reference.somalier_sites_vcf,
        runtime_attributes = default_runtime_attributes
    }
  }

  # check for sample swaps and confirm gender of sample
  call Somalier.somalier_relate {
    input:
      sample_id = sample.sample_id,
      extracted_somalier_sites = somalier_extract.extracted_sites,
      runtime_attributes = default_runtime_attributes
  }

  if (somalier_relate.min_relatedness >= min_relatedness_sample_swap) {
    
    if (somalier_relate.inferred_sex != "UNKNOWN") {

      scatter (idx in range(length(pbmm2_align.bam))) {
        # call structural variants with sniffles
        call Sniffles.sniffles_discover {
          input:
            sample_id = sample.sample_id,
            bam = aligned_bam[idx].data,
            bam_index = aligned_bam[idx].index,
            reference_fasta = reference.fasta.data,
            reference_index = reference.fasta.index,
            tr_bed = reference.tandem_repeat_bed,
            runtime_attributes = default_runtime_attributes
        }
      }

      scatter (idx in range(length(aligned_bam))) {
        # discover SV signatures with pbsv
        call Pbsv.pbsv_discover {
          input:
            bam = aligned_bam[idx].data,
            bam_index = aligned_bam[idx].index,
            tr_bed = reference.tandem_repeat_bed,
            runtime_attributes = default_runtime_attributes
        }
      }
      Array[File] svsig = flatten(pbsv_discover.svsigs)

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

      # merge aligned bams for trgt
      call Utils.merge_bams {
        input:
          bams = pbmm2_align.bam,
          output_bam_name = "~{sample.sample_id}.~{reference.name}.bam",
          runtime_attributes = default_runtime_attributes
      }

      IndexData merged_bam = {
        "data": merge_bams.merged_bam, 
        "index": merge_bams.merged_bam_index
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

      # keep track of sample IDs that pass QC
      String qc_pass = sample.sample_id
    }
  }

  if (somalier_relate.min_relatedness < min_relatedness_sample_swap) {
    # keep track of sample IDs that are excluded for sample swap between movies
    String qc_fail_relatedness = sample.sample_id
  }
  if (somalier_relate.inferred_sex == "UNKNOWN") {
    # keep track of sample IDs that are excluded for unknown sex
    String qc_fail_sex = sample.sample_id
  }

  output {
    # QC pass
    IndexData? pbmm2_merged_bam = merged_bam
    Array[File]? sniffles_snfs = sniffles_discover.snf
    Array[File]? pbsv_svsigs = svsig
    IndexData? deepvariant_gvcf = deepvariant.gvcf
    IndexData? trgt_vcf = { "data": select_first([trgt.repeat_vcf]), "index": select_first([trgt.repeat_vcf_index]) }
    String? qc_pass_id = qc_pass
    
    # QC fail
    String? qc_fail_relatedness_id = qc_fail_relatedness
    String? qc_fail_sex_id = qc_fail_sex
  }
}
