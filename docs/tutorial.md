# Running the workflow

**Please [contact Juniper Lake by email](mailto:jlake@pacificbiosciences.com) before attempting to run this workflow. This is NOT a development-level workflow and requires special instructions.**

## HPC Quickstart

This sections describes running the workflow with MiniWDL on a SLURM job scheduler.

### Requirements

- Your HPC job scheduler should be SLURM
- Python3
- Pip

> [!WARNING]
> Each step should be completed without errors before moving on to the next.

### 1. Download and install everything

```
# clone github repo
git clone https://github.com/juniper-lake/CoLoRSdb.git

# make virtual environment and install dependencies
python3 -m venv .venv
./.venv/bin/pip install -U pip
./.venv/bin/pip install miniwdl>=1.9.1 miniwdl-slurm

# download and unzip required reference files
wget https://zenodo.org/records/10277930/files/colorsdb.v1.0.1.resources.tgz
```

### 2. Create your input/configuration files

There are two files that need to be copied to the working directory and edited to specify the correct inputs and workflow configuration.

First, copy them to your working directory. Don't move these files from their original locations because some are used for testing.

```
cp CoLoRSdb/wdl/tests/test_data/sample_sheet.tsv CoLoRSdb/docs/sources/miniwdl.cfg .
```

Second, edit the files with your favorite editor.

- Update `miniwdl.cfg` so it will run correctly on your system. Please see the [default miniwdl config](https://github.com/chanzuckerberg/miniwdl/blob/main/WDL/runtime/config_templates/default.cfg) and the [miniwdl-slurm config example](https://github.com/miniwdl-ext/miniwdl-slurm#configuration) for more details.
- Replace sample info in `sample_sheet.tsv` with your own. First column is sample IDs, which should have no spaces or special characters except underscores. The second column is a comma-separated list of HiFi movies (FASTQ or BAM) associated with the sample.

### 3. Test to make sure miniwdl works

> [!WARNING]
> You cannot launch miniwdl on SLURM from an interactive `srun` session.

```
# activate your virtual environment
source .venv/bin/activate

# test miniwdl and miniwdl-slurm, this should complete quickly
miniwdl run --verbose \
  --dir miniwdl_execution/tests \
  --cfg miniwdl.cfg \
  CoLoRSdb/wdl/tests/hello.wdl
```

### 4. Run the workflow

```
# activate your virtual environment if not already activated
source .venv/bin/activate

# run your workflow on GRCh38
# replace anything in <> with your own info
miniwdl run --verbose \
  --cfg miniwdl.cfg \
  --dir miniwdl_execution \
  --input inputs.hpc.grch38.json \
  CoLoRSdb/wdl/workflows/colors_main.wdl \
  cohort_id=<your_cohort_id> \
  sample_sheet=<path/to/sample_sheet.tsv> \
  reference_bundle=colorsdb.v1.0.1.resources.tgz \
  backend=HPC \
  preemptible=false
```

## AnViL/Terra Quickstart

Coming soon!

  <!-- "colors_main.backend": "AnVIL",
  "colors_main.preemptible": true,
  "colors_main.zones": "us-central1-a us-central1-c us-central1-b us-central1-f" -->
