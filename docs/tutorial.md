# Running the workflow

**Workflow entrypoint**: [wdl/workflows/main.wdl](wdl/workflows/main.wdl)

While this workflow is structured to be run a variety of backends (HPC, AWS, GCP, Azure) using either `Miniwdl` or `Cromwell`, it has only been tested on HPC (SLURM + MiniWDL) and AnViL (GCP/Azure + Cromwell). The following instructions are limited to these two options.

- [HPC Quickstart](#hpc-quickstart)
- [AnViL/Terra Quickstart](#anvilterra-quickstart)

## HPC Quickstart

This sections describes running the workflow with MiniWDL on a SLURM job scheduler.

### Requirements

- Your HPC job scheduler should be SLURM
- Python3
- Pip

> [!WARNING] Each step should be completed without errors before moving on to the next.

### 1. Download and install everything

```
# clone github repo
git clone https://github.com/juniper-lake/CoLoRSdb.git

# make virtual environment and install dependencies
cd CoLoRSdb
python3 -m venv .venv
./.venv/bin/pip install -U pip
./.venv/bin/pip install miniwdl>=1.9.1 miniwdl-slurm

# download and unzip required reference files
wget https://zenodo.org/records/10277930/files/colorsdb.v1.0.1.resources.tgz
tar -xzf colorsdb.v1.0.1.resources.tgz && rm colorsdb.v1.0.1.resources.tgz
```

### 2. Create your input/configuration files

There are four files that need to be copied to the working directory and edited to specify the correct inputs and workflow configuration.

First, copy them to your working directory. Don't move these files from their original locations because some are used for testing.

```
cp wdl/tests/test_data/sample_sheet.tsv backends/hpc/* .
```

Second, edit the files with your favorite editor.

- Update `miniwdl.cfg` so it will run correctly on your system. Please see the [default miniwdl config](https://github.com/chanzuckerberg/miniwdl/blob/main/WDL/runtime/config_templates/default.cfg) and the [miniwdl-slurm config example](https://github.com/miniwdl-ext/miniwdl-slurm#configuration) for more details.
- Update both `inputs.hpc.grch38.json` and `inputs.hpc.chm13.json` with the correct `cohort_id` and `sample_sheet`. These values should be the same in both files.
- Replace sample info in `sample_sheet.tsv` with your own. First column is sample IDs, which should have no spaces or special characters except underscores. The second column is a comman-separated list of HiFi movies (FASTQ or BAM) associated with the sample.

### 3. Test to make sure miniwdl works

> [!WARNING]
> You cannot launch miniwdl from an interactive `srun` session.

```
# activate your virtual environment
source .venv/bin/activate

# test miniwdl and miniwdl-slurm, this should complete quickly
miniwdl run --verbose \
  --dir miniwdl_execution/tests \
  --cfg miniwdl.cfg \
  wdl/tests/hello.wdl
```

### 4. Run the workflow

Remember to run on both GRCh38 and CHM13.

```
# activate your virtual environment if not already activated
source .venv/bin/activate

# run your workflow on GRCh38
miniwdl run --verbose \
  --cfg miniwdl.cfg \
  --dir miniwdl_execution/grch38 \
  --input inputs.hpc.grch38.json \
  wdl/workflows/main.wdl \

# when previous run is complete, run on CHM13
miniwdl run --verbose \
  --cfg miniwdl.cfg \
  --dir miniwdl_execution/chm13 \
  --input inputs.hpc.chm13.json \
  wdl/workflows/main.wdl \
```

## AnViL/Terra Quickstart

Coming soon!
