# Workflow outputs

## Cohort-level outputs

These files will be output for each cohort if the listed requirements are met.

\* QC = For these outputs to be produced, at least 2 samples must pass quality control.

### QC outputs

| Type | Name | Description | Requirements |
| :- | :- | :- | :- |
| String | sample_size_message | Message conveying insufficient sample size when the number of samples is less than 2 | |
| String | qc_message | Message conveying insufficient sample size when the number of samples passing QC is less than 2 | |
| File? | qc_summary_tsv | TSV of various quality metrics for each sample | Number of samples > 1  |
| File? | somalier_pairs_tsv | TSV of pairwise relatedness, "pairs.tsv" from [somalier](https://github.com/brentp/somalier) | QC* |
| File? | somalier_samples_tsv | "samples.tsv" from [somalier](https://github.com/brentp/somalier) | QC* |

### Zipped, indexed, multi-sample VCFs + TBI index

If `anonymize_output = false`, both the original and postprocessed VCF (with ploidy corrected on sex chromosomes) are included as outputs.

| Type | Name | Description | Requirements |
| :- | :- | :- | :- |
| [IndexData](wdl/structs.wdl)? | deepvariant_glnexus_postprocessed_vcf | Small variants (SNVs, small indels) called by [DeepVariant](https://github.com/google/deepvariant) and joint called with [GLnexus](https://github.com/dnanexus-rnd/GLnexus) | QC* |
| Array[[IndexData](wdl/structs.wdl)]? | pbsv_jasminesv_postprocessed_vcf | Structural variants (only INS,DEL,INV types) called with [pbsv](https://github.com/PacificBiosciences/pbsv) and merged with [jasmineSV](https://github.com/mkirsche/Jasmine). | QC* |
| [IndexData](wdl/structs.wdl)? | sniffles_postprocessed_vcf | Structural variants joint called with [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) | QC* |
| [IndexData](wdl/structs.wdl)? | trgt_postprocessed_vcf | Tandem repeats genotyped with [TRGT](https://github.com/PacificBiosciences/trgt) and merged using a [custom script](docker/vcfparser/scripts/merge_trgt_vcfs.py) | QC* ; `reference.trgt_tandem_repeat_bed` |
| [IndexData](wdl/structs.wdl)? | hificnv_postprocessed_vcf | CNVs called with [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) and merged with [bcftools](https://github.com/samtools/bcftools) | QC* ; `reference.hificnv_exclude_bed` ; `reference.hificnv_expected_bed_male` ; `reference.hificnv_expected_bed_female` |

### VCF stats

| Type | Name | Description | Requirements |
| :- | :- | :- | :- |
| Array[File]? | pbsv_vcf_stats | Per-sample variant statistics from `cohort_pbsv_vcf` queried with [bcftools](https://github.com/samtools/bcftools) | QC* |
| File? | sniffles_vcf_stats | Per-sample variant statistics from `cohort_sniffles_vcf` queried with [bcftools](https://github.com/samtools/bcftools) | QC* |

### Ancestry outputs

| Type | Name | Description | Requirements |
| :- | :- | :- | :- |
| File? | peddy_het_check | A CSV of samples with higher levels of HET calls, lower depth, or more variance in b-allele-frequency (ref / (ref + alt )) for het calls, produced by [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_sex_check | A CSV of discrepancies between ped-reported and genotype-inferred sex, produced by [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_background_pca | An ancestry prediction based on projection onto the thousand genomes principal components produced by [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_html | Interactive html output from [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_ped | A ped file that also lists the most useful columns from the het-check and sex-check, produced by [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
| File? | peddy_vs_html | Interactive html output from [peddy](https://github.com/brentp/peddy) | QC* ; `reference.peddy_sites` ; `reference.peddy_bin` |
