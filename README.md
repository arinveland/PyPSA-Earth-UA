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

## Choosing a solver
The following open source solvers are pre-installed in the PyPSA-Earth environment: GLPK, WinGLPK, HiGHS. Gurobi is also pre-installed, but you must provide your own licence to use it.