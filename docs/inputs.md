# Workflow inputs

Reference datasets and associated input files are hosted publicly for use in the pipeline. 

You can validate your input json file against the [schema](../wdl/workflows/main.input.schema.json) using [`check-jsonschema`](https://check-jsonschema.readthedocs.io/en/latest/index.html) (requires `python3` and `pip`). This will ensure that inputs are the correct format and that no required inputs are missing.

```bash
# install required packages
pip install check-jsonschema 'jsonschema[format]'

# validate json
check-jsonschema --schemafile ./wdl/workflows/main.input.schema.json <your_inputs.json> 
```

## Cohort-specific inputs

A cohort should include at least 2 samples, although if you're running this on only 2 samples, this probably isn't the workflow for you. Samples can either be specified using **either** `sample_sheet` (example [here](examples/sample_sheet.tsv)) **or** by defining both `sample_ids` and `movies`. 

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | cohort_id | A unique name for the cohort; used to name outputs | |
| Array[String]+? | sample_ids | The set of samples for the cohort. When paired with `movies`, this can be used as an alternative input to `sample_sheet`. | Preferred sample input for Terra/AnViL. |
| Array[Array[File]]+? | movies | The movies for each sample. When paired with `sample_ids`, this can be used as an alternative input to `sample_sheet` | Preferred sample input for Terra/AnViL. |
| File? | sample_sheet | A sample sheet specifying `sample_ids` and their associated `movies` for the cohort. See an [example](examples/sample_sheet.tsv) | Preferred sample input for HPC. |

## Analysis inputs

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| Boolean | anonymize_output | Shuffle genotypes in output VCFs to de-identify samples | \[true, false\] |
| Int | max_samples_pbsv_call | Maximum number of samples to jointly call with pbsv, used to chunk the pbsv_call task | \[2,inf\) |
| Float | max_sample_relatedness_qc | Maximum relatedness between any two samples; samples will be iteratively removed from the cohort to satisfy this quality control threshold | \[0,1\] |
| Float | min_movie_relatedness_qc | Minimum relatedness between any two movies from the same sample; samples where all movies do not satisfy this quality conrol threshold will be removed from the cohort | \[0,1\]

## [ReferenceData](workflows/humanwgs_structs.wdl)

Inputs associated with the reference genome. All files are hosted publicly on zenodo. Only **GRCh38** and **CHM13** reference genomes are currently supported for this workflow, and some workflow steps (TRGT, HiFiCNV, peddy) only support GRCh38 (see notes below).

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | name | Reference name; used to name outputs (e.g., "GRCh38") | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | fasta | Reference genome and index | |
| Array[String] | chromosomes | Chromosomes in which to call structural variants with pbsv (for parallelization) | |
| File | non_diploid_regions | BED-style list of regions where ploidy < 2 with added columns for sex and ploidy, used to correct genotype ploidy in output VCFs | |
| File | tandem_repeat_bed | Tandem repeat locations used by [pbsv](https://github.com/PacificBiosciences/pbsv) and [sniffles](https://github.com/fritzsedlazeck/Sniffles) to normalize SV representation | |
| File | somalier_sites_vcf | Sites used for quality control on aligned BAMs with [somalier](https://github.com/brentp/somalier) | |
| File? | trgt_tandem_repeat_bed | Tandem repeat sites to be genotyped by [TRGT](https://github.com/PacificBiosciences/trgt) | GRCh38 only ; required for [TRGT](https://github.com/PacificBiosciences/trgt) to run |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | hificnv_exclude_bed | Compressed BED and index of regions to exclude from calling by [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV). | GRCh38 only ; required for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) to run |
| File? | hificnv_expected_bed_male | BED of expected copy number for male karyotype for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) | GRCh38 only ; required for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) to run |
| File? | hificnv_expected_bed_female | BED of expected copy number for female karyotype for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) | GRCh38 only ; required for [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) to run |
| File? | peddy_sites | sites to be used for ancestry estimation with [peddy](https://github.com/brentp/peddy) | GRCh38 only ; required for [peddy](https://github.com/brentp/peddy) to run |
| File? | peddy_bin |  raw binary alternate counts (gt_types) from thousand-genomes that have been written as uint8 and gzipped for ancestry estimation with [peddy](https://github.com/brentp/peddy) | GRCh38 only ; required for [peddy](https://github.com/brentp/peddy) to run |

## Runtime inputs

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | deepvariant_version | Override the default DeepVariant software version | |
| Int? | pbsv_call_mem_gb | Optionally set RAM (GB) for pbsv_call joint calling | |
| Int? | glnexus_mem_gb | Optionally set RAM (GB) for GLnexus joint calling | |
| Int? | sniffles_call_mem_gb | Optionally set RAM (GB) for sniffles joint calling | |
| String | backend | Backend where the workflow will be executed | \["Azure", "AWS", "GCP", "HPC"\] |
| String? | zones | Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'. | <ul><li>[Determining available zones in AWS](backends/aws/README.md#determining-available-zones)</li><li>[Determining available zones in GCP](backends/gcp/README.md#determining-available-zones)</li></ul> |
| String? | aws_spot_queue_arn | Queue ARN for the spot batch queue; required if backend is set to 'AWS' and `preemptible` is set to `true` | [Determining the AWS queue ARN](backends/aws/README.md#determining-the-aws-batch-queue-arn) |
| String? | aws_on_demand_queue_arn | Queue ARN for the on demand batch queue; required if backend is set to 'AWS' and `preemptible` is set to `false` | [Determining the AWS queue ARN](backends/aws/README.md#determining-the-aws-batch-queue-arn) |
| Boolean | preemptible | If set to `true`, run tasks preemptibly where possible. On-demand VMs will be used only for tasks that run for >24 hours if the backend is set to GCP. If set to `false`, on-demand VMs will be used for every task. Ignored if backend is set to HPC. | \[true, false\] |
| String | container_registry | The repo where docker images are hosted | |

