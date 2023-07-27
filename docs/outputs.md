# Workflow outputs

## Cohort-level outputs

These files will be output for each cohort if the listed requirements are met.

| Type | Name | Description | Requirements |
| :- | :- | :- | :- |
| String | message | Message conveying insufficient sample size when the number of samples is less than 2 | |
| File? | qc_summary_tsv | TSV of various quality metrics for each sample | Number of samples > 1  |
| File? | pairwise_relatedness_tsv | TSV of pairwise relatedness, "pairs.tsv" from [somalier](https://github.com/brentp/somalier) | QC* |
| [IndexData](wdl/structs.wdl)? | cohort_deepvariant_vcf | Multi-sample VCF+index of small variants (SNVs, small indels) called by [DeepVariant](https://github.com/google/deepvariant) and joint called with [GLnexus](https://github.com/dnanexus-rnd/GLnexus) | QC* |
| Array[[IndexData](wdl/structs.wdl)]? | cohort_pbsv_vcf | Multi-sample VCFs+indexes of structural variants (only INS,DEL,INV types) joint called with [pbsv](https://github.com/PacificBiosciences/pbsv). Multiple VCFs are produced if the number of samples that pass QC exceeds the input `max_samples_pbsv_call` | QC* |
| [IndexData](wdl/structs.wdl)? | cohort_sniffles_vcf | Multi-sample VCF+index of structural variants joint called with [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) | QC* |
| [IndexData](wdl/structs.wdl)? | cohort_trgt_vcf | Multi-sample VCF+index of tandem repeats genotyped with [TRGT](https://github.com/PacificBiosciences/trgt) and merged using a [custom script](docker/vcfparser/scripts/merge_trgt_vcfs.py) | QC* ; `reference.trgt_tandem_repeat_bed` |
| [IndexData](wdl/structs.wdl)? | cohort_hificnv_vcf | Multi-sample VCF+index of CNVs called with [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) and merged with [bcftools](https://github.com/samtools/bcftools) | QC* ; `reference.hificnv_exclude_bed` ; `reference.hificnv_expected_bed_male` ; `reference.hificnv_expected_bed_female` |
| File? | cohort_deepvariant_vcf_stats | Per-sample variant statistics from `cohort_deepvariant_vcf` produced with [bcftools stats](https://github.com/samtools/bcftools) | QC* |
| Array[File]? | cohort_pbsv_vcf_stats | Per-sample variant statistics from `cohort_pbsv_vcf` queried with [bcftools](https://github.com/samtools/bcftools) | QC* |
| File? | cohort_sniffles_vcf_stats | Per-sample variant statistics from `cohort_sniffles_vcf` queried with [bcftools](https://github.com/samtools/bcftools) | QC* |
| File? | peddy_het_check | A CSV of samples with higher levels of HET calls, lower depth, or more variance in b-allele-frequency (ref / (ref + alt )) for het calls, produced by [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_sex_check | A CSV of discrepancies between ped-reported and genotype-inferred sex, produced by [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_background_pca | An ancestry prediction based on projection onto the thousand genomes principal components produced by [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_html | Interactive html output from [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_ped | A ped file that also lists the most useful columns from the het-check and sex-check, produced by [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_vs_html | Interactive html output from [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |

\* QC = For these outputs to be produced, at least 2 samples must pass quality control.