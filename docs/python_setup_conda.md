# Setting up a Python environment using `conda`

You can do this on JASMIN, or on your local computer. If on JASMIN,
log in to the notebook service: https://notebooks.jasmin.ac.uk/

Open a terminal in JupyterLab: File → New → Terminal

Clone this repo:

```
git clone https://github.com/digital-earths-UK-hackathon/hk26.git
cd hk26
```

Create a Python env and add the kernel to the notebook service:

```
conda env create -f conda_envs/hk26_env.yaml
conda activate hk26_env
python -m ipykernel install --user --name hk26_env --display-name "Hackathon 2026 (conda)"
```

You can now create a new notebook with the "Hackathon 2026 (conda)" kernel to use all of the packages.

To install additional packages:

```
conda activate hk26_env
conda install <package_name>
```

## JASMIN help docs

See [here for more information about using the JASMIN notebook service](https://help.jasmin.ac.uk/docs/interactive-computing/jasmin-notebooks-service/).

## Running on JASMIN ssh

If you are `ssh`ing into JASMIN and want to run on the command line, you may need to run `module load jaspy` first before
using and `conda` commands or activating your environment.