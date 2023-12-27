# Workflow inputs

## Cohort-specific inputs

A cohort should include at least 2 samples, although if you're running this on only 2 samples, this probably isn't the workflow for you.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | cohort_id | A unique name for the cohort; used to name outputs | |
| File | sample_sheet | A sample sheet specifying `sample_ids` and their associated `movies` for the cohort. | See an [example](examples/sample_sheet.tsv) |
| Boolean | anonymize_output | Shuffle genotypes in output VCFs to de-identify samples | \[true, false\] |

## [ReferenceData](workflows/humanwgs_structs.wdl)

Inputs associated with the reference genome. All files are hosted publicly on zenodo. Only **GRCh38** and **CHM13** reference genomes are currently supported for this workflow, and some workflow steps (TRGT, HiFiCNV, peddy) only support GRCh38 (see notes below).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10277930.svg)](https://doi.org/10.5281/zenodo.10277930)

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | name | Reference name; used to name outputs (e.g., "GRCh38") | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | fasta | Reference genome and index | |
| Array[String] | chromosomes | Chromosomes in which to call structural variants with pbsv (for parallelization) | |
| File | non_diploid_regions | BED-style list of regions where ploidy < 2 with added columns for sex and ploidy, used to correct genotype ploidy in output VCFs | |
| File | tandem_repeat_bed | Tandem repeat locations used by [pbsv](https://github.com/PacificBiosciences/pbsv) and [sniffles](https://github.com/fritzsedlazeck/Sniffles) to normalize SV representation | |
| File | somalier_sites_vcf | Sites used for quality control on aligned BAMs with [somalier](https://github.com/brentp/somalier) | |
| Array[File]? | trgt_tandem_repeat_beds | Tandem repeat sites to be genotyped by [TRGT](https://github.com/PacificBiosciences/trgt) | GRCh38 only ; required for [TRGT](https://github.com/PacificBiosciences/trgt) to run |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | hificnv_exclude_bed | Compressed BED and index of regions to exclude from calling by [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV). | GRCh38 only ; required for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) to run |
| File? | hificnv_expected_bed_male | BED of expected copy number for male karyotype for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) | GRCh38 only ; required for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) to run |
| File? | hificnv_expected_bed_female | BED of expected copy number for female karyotype for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) | GRCh38 only ; required for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) to run |
| File? | peddy_sites | sites to be used for ancestry estimation with [peddy](https://github.com/brentp/peddy) | GRCh38 only ; required for [peddy](https://github.com/brentp/peddy) to run |
| File? | peddy_bin |  raw binary alternate counts (gt_types) from thousand-genomes that have been written as uint8 and gzipped for ancestry estimation with [peddy](https://github.com/brentp/peddy) | GRCh38 only ; required for [peddy](https://github.com/brentp/peddy) to run |

## Runtime inputs

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | backend | Backend where the workflow will be executed | \["AnVIL", "Azure", "AWS", "GCP", "HPC"\] |
| String? | zones | Zones where compute will take place; required if backend is set to 'AWS', 'GCP', or 'AnVIL'. | |
| String? | aws_spot_queue_arn | Queue ARN for the spot batch queue; required if backend is set to 'AWS' and `preemptible` is set to `true` | |
| String? | aws_on_demand_queue_arn | Queue ARN for the on demand batch queue; required if backend is set to 'AWS' and `preemptible` is set to `false` | |
| Boolean | preemptible | If set to `true`, run tasks preemptibly where possible. On-demand VMs will be used only for tasks that run for >24 hours if the backend is set to GCP. If set to `false`, on-demand VMs will be used for every task. Ignored if backend is set to HPC. | \[true, false\] |
