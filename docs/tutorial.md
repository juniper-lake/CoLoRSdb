# Running the workflow

**Workflow entrypoint**: [wdl/workflows/main.wdl](wdl/workflows/main.wdl)

## Setting up and running the workflow

1. Install and configure the workflow execution engine of your choice following the documentation for the backend environment where your data is located. See the [backend environments](#backend-environments) section for more information specific to your backend.

2. Fill out the [input template file](../wdl/workflows/input_template.json). Partially pre-filled input templates can also be found in [backends](../backends).

3. Optionally validate your input json by following the instructions in [inputs](inputs.md).

4. Run the workflow using the engine and backend of choice (described below). Please read backend-specific information on launching the workflow in [backends](../backends).

## Backend environments

The workflow can be run on Azure, AWS, GCP, or HPC, but has thus far only been tested on HPC and AnVIL (GCP). 

- [AnVIL](backends/anvil)
- [HPC](backends/hpc)

## Workflow engines

Two popular engines for running WDL-based workflows are [`miniwdl`](https://miniwdl.readthedocs.io/en/latest/getting_started.html) and [`Cromwell`](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

The workflow engine that you choose will depend on where your data is located.

| Engine | Azure | AWS | GCP | HPC |
| :- | :- | :- | :- | :- |
| [**miniwdl**](https://github.com/chanzuckerberg/miniwdl#scaling-up) | _Unsupported_ | Supported via the [Amazon Genomics CLI](https://aws.amazon.com/genomics-cli/) | _Unsupported_ | (SLURM only) Supported via the [`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm) plugin |
| [**Cromwell**](https://cromwell.readthedocs.io/en/stable/backends/Backends/) | Supported via [Cromwell on Azure](https://github.com/microsoft/CromwellOnAzure) | Supported via the [Amazon Genomics CLI](https://aws.amazon.com/genomics-cli/) | Supported via Google's [Pipelines API](https://cromwell.readthedocs.io/en/stable/backends/Google/) | Supported - [Configuration varies depending on HPC infrastructure](https://cromwell.readthedocs.io/en/stable/tutorials/HPCIntro/) |

### Run using miniwdl

`miniwdl run wdl/workflows/main.wdl --input <input_file_path.json>`

### Run using Cromwell

`java -jar <cromwell_jar_path> run wdl/workflows/main_cohort.wdl -i <input_file_path.json>`
