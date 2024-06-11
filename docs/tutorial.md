# Running the workflow

**Please [contact Juniper Lake by email](mailto:jlake@pacificbiosciences.com) before attempting to run this workflow or if you hit any errors while running the workflow. If you are contributing data to CoLoRSdb, then customizations to the workflow must be approved to ensure the analysis is the same for each contributor.**

## HPC Quickstart (SLURM+MiniWDL or Cromwell)

If running with SLURM+MiniWDL, the following requirements must be met.

- Your HPC job scheduler should be SLURM
- Python3
- Pip

> [!WARNING]
> Each step should be completed without errors before moving on to the next.

### 1. Download and install everything

```
# clone github repo
git clone https://github.com/juniper-lake/CoLoRSdb.git

# download required reference files (do not unzip)
wget https://zenodo.org/records/10277930/files/colorsdb.v1.0.1.resources.tgz

# if using miniwdl, make virtual environment and install dependencies
python3 -m venv .venv
./.venv/bin/pip install -U pip
./.venv/bin/pip install miniwdl>=1.9.1 miniwdl-slurm
```

### 2. Create your input/configuration files

Create your own [sample sheet TSV](templates/sample_sheet.tsv) and [inputs JSON](templates/input_template.json) based on the linked templates. Detailed description of inputs are below.

| Input | Type | Description |
| --- | --- | --- |
| cohort_id | String | Name of the cohort, e.g. "HPRC" |
| sample_sheet | File | TSV where first column is sample IDs, which should have no spaces or special characters except underscores. The second column is a comma-separated list of HiFi movies (FASTQ or BAM) associated with the sample. BAMs can be aligned or unaligned |
| reference_bundle | File | Zipped tarball of reference files downloaded according to above instructions, i.e. "colorsdb.v1.0.1.resources.tgz" |
| anonymize_output | Boolean | Default is `false`. Set to `true` if you CANNOT share sample-level data with PacBio. If `true` then we cannot use your variant data for sex chromosomes because only randomized data will be output by the workflow. See how data is randomized [here](images/anonymize_output_example.png) |
| backend | String | Can be "HPC", "AnViL", "AWS", "Azure", or "GCP" depending on your backend |
| preemptible | Boolean | Set to `true` to run tasks preemptibly where possible. Set to `false` to use on-demand VMs for every task. Ignored if backend is set to HPC |

### 3. Configure and test MiniWDL (Cromwell users can skip)

- Create a [miniwdl config file](templates/miniwdl.cfg) using the linked template. For additional miniwdl configuration details, please see the [default miniwdl config](https://github.com/chanzuckerberg/miniwdl/blob/main/WDL/runtime/config_templates/default.cfg) and a [miniwdl-slurm config example](https://github.com/miniwdl-ext/miniwdl-slurm#configuration).

```
# activate your virtual environment
source .venv/bin/activate

# test miniwdl and miniwdl-slurm, this should complete quickly
# replace anything in <> with your own info
miniwdl run --verbose \
  --dir miniwdl_execution/tests \
  --cfg <path/to/your/miniwdl.cfg> \
  CoLoRSdb/wdl/tests/hello.wdl
```

### 4. Run the workflow

#### MiniWDL

> [!WARNING]
> You cannot launch miniwdl on SLURM from an interactive `srun` session or from an `sbatch` job. You must run from a login node.

```
# activate your virtual environment if not already activated
source .venv/bin/activate

# replace anything in <> with your own info
miniwdl run --verbose \
  --cfg <path/to/your/miniwdl.cfg> \
  --dir miniwdl_execution \
  --input <path/to/your/inputs.json> \
  CoLoRSdb/wdl/workflows/colors_main.wdl
```

#### Cromwell

```
# replace anything in <> with your own info
cromwell run \
  CoLoRSdb/wdl/workflows/colors_main.wdl \
  --inputs <path/to/your/inputs.json>
```
