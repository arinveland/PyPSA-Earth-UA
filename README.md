# Running the model

## Installing Dependencies
It is recommended to install dependencies using a package manager such as Conda. For supported processors and operating systems, it is recommended to use the corresponding locked files. In other cases the unlocked environment file may be used, though this may come with compatibility issues.

### Installation and activation with Conda
Run the following commands in your terminal:
```console
conda env create -f Files/envs/<INSERT FILENAME>.lock.yaml
```

```console
conda activate pypsa-earth
```

### Installation and activation with Mamba
Run the following commands in your terminal:
```console
conda install -c conda-forge mamba
```

```console
mamba env create -f Files/envs/<INSERT FILENAME>.lock.yaml
```

```console
conda activate pypsa-earth
```

### Deactivation
To deactivate the environment, run the following command in your terminal:
```console
conda deactivate
```

## Choose a solver
The following open source solvers are pre-installed in the PyPSA-Earth environment: GLPK, WinGLPK, HiGHS. Gurobi is also pre-installed, but you must provide your own licence to use it. By default, the scripts in this workflow try Gurobi and HiGHS first, as these are the most capable ones. This can be changed by altering the CANDIDATE_SOLVERS list in the scripts.

## Add infrastructure attack data
This should be added on csv format to the Files folder. As a suggestion, ACLED has a dataset sourced from open media reports that they are often willing to share with serious researchers (https://acleddata.com/monitor/ukraine-conflict-monitor).

## Run the model
To run any of the snakemake jobs in the workflow, use the command:

```console
snakemake -j [cores] [job_name]
```

Where cores is the number of provided cores, and job_name is the name of the snakemake job you want to run (ua_scenario_analysis to run a scenario). Other jobs will then be run as required based on which inputs already exist.

This model is quite computationally heavy, espeacially when run over a long time horizon, so it is recommended to run it on a HPC cluster.