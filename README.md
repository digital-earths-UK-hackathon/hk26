# Repository for the UK km-scale hackathon 2026

## Getting started on JASMIN

If you haven't already set up Python on JASMIN, see [this guide](docs/python_setup_conda.md).

If planning to commit any code to community repositories on GitHub, you may want to check [this GitHub setup guide](docs/github_setup.md)

## Dataset catalog

You can then open the datasets contained in [this intake 0.7 catalog](https://digital-earths-global-hackathon.github.io/catalog/).
You can find all Unified Model by looking for entried that start with `um_`, and those that were created for this hackathon
by their `_hk26` suffix.

## K-Scale R&D datasets availability and access

For the DYAMOND-3 year long simulations the native .pp outputs from the Met Office Unified Model are currently available on JASMIN via the kscale group workspace. If you are not already part of the group workspace, please follow these [instructions](https://help.jasmin.ac.uk/docs/short-term-project-storage/apply-for-access-to-a-gws/) and search kscale for the group workspace name. Below are the paths to the datasets:

**Global 5km RAL3.3** (explicit convection [see Bush et al., 2025](https://gmd.copernicus.org/articles/18/3819/2025/)): /gws/ssde/j25b/kscale/DYAMOND3_reruns/5km-RAL3p3-tuned/glm/field.pp/
**Global 5km CoMA9_TBv1.2** (CoMorph-A convection scheme with changes for running in the convective grey-zone): /gws/ssde/j25b/kscale/DYAMOND3_reruns/5km-CoMA9/glm/field.pp/
**Global 10km CoMA9_TBv1.2** (science configuration as above): /gws/ssde/j25b/kscale/DYAMOND3_reruns/10km-CoMA9/glm/field.pp/
**Global 10km GAL9** (6A convection scheme see [Willett et al., 2025](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-1829/)): /gws/ssde/j25b/kscale/DYAMOND3_reruns/10km-GAL9/glm/field.pp/

As part of the K-Scale model hierarchy we also run a cyclic tropical channel (CTC) and continental scale limited area models (LAMs). The currently available CTC and LAM datasets are nested inside the Global 5km CoMA9_TBv1.2 simulation. For the CTC and each of the three LAM domains both the RAL3.3 and CoMA9_TBv1.2 science configurations are tested, with the horizontal grid spacing of these simulations 0.04 degrees (~4.4 km). In the short term the CTC and LAM datasets are stored on scratch on JASMIN, this page will be updated when the datasets are transferred to MASS storage. 

CTC 4.4 km CoMA9_TBv1.2: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/CTC_km4p4_CoMA9_TBv1/field.pp/
CTC 4.4 km RAL3.3: Not yet available on JASMIN

South America LAM 4.4 km CoMA9_TBv1.2: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/SAmer_km4p4_CoMA9_TBv1/field.pp/
South America LAM 4.4 km RAL3.3: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/SAmer_km4p4_RAL3P3/field.pp/
Africa LAM 4.4 km CoMA9_TBv1.2: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/Africa_km4p4_CoMA9_TBv1/field.pp/
Africa LAM 4.4 km RAL3.3: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/Africa_km4p4_RAL3P3/field.pp/
SE Asia LAM 4.4 km CoMA9_TBv1.2: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/SEA_km4p4_CoMA9_TBv1/field.pp/
SE Asia LAM 4.4 km RAL3.3: /work/scratch-pw5/rwjones/kscale/DYAMOND3_reruns/5km-CoMA9/SEA_km4p4_RAL3P3/field.pp/

The Global 5km demonstrator simulations (and Global 10km control) are 15-day forecasts re-initialised every 3 days for the period 1 September 2020 to end August 2021. The pp output from these simulations are stored on MASS and can be accessed via JASMIN - the relevant suite IDs for accessing the dataset are below:

Global 5km CoMA9_TBv1.2: moose:/devfc/u-dp135/field.pp/ - with files of the format YYYYMMDDT0000Z_umglaa_p*.pp
Global 10km CoMA9_TBv1.2: moose:/devfc/u-dq613/field.pp/ - with files of the format YYYYMMDDT0000Z_umglaa_p*.pp
Global 10km GAL9 control: suite ID to be added soon

This page will soon be updated with details of how to access the 40-day DYAMOND-Summer and DYAMOND-Winter Unified Model experiments (via MASS). 

## JASMIN notebooks

There are some notebooks that give examples of how to access the data, for example [accessing global data](notebooks/01_view_global_data.ipynb).
To use these on JASMIN, clone this repo into your home directory and then navigate to the `notebooks` directory in the
[JASMIN notebook service](https://notebooks.jasmin.ac.uk/).
Or you can directly download the notebook and upload to the [JASMIN notebook service](https://notebooks.jasmin.ac.uk/).
The notebooks can also be run locally; just set use the `online` catalog.

## Sharing code

There are directories for the different groups in this repo, all prefixed with `hk26-`. Any code and notebooks can be shared here, for notebooks please
check they run with the standard conda env and clear all output before committing. 
You cannot push directly to the main branch, instead, follow these instructions to create a Pull Request PR:

```bash
# Make sure main is up-to-date
git checkout main
git pull

# Create a new branch locally with a sensible branch name
git checkout -b my_new_branch

# Commit your changes
git status
git add .
git commit

# Push your changes to github
# If you get an error, email your github username to Mark Muetzelfeldt, email below.
git push origin my_new_branch

# Go to https://github.com/digital-earths-UK-hackathon/hk26 and look for the option to create new pull request.
# Or go straight to https://github.com/digital-earths-UK-hackathon/hk26/pulls
# Create and merge your PR into main.
# You can add reviewers here - this is recommended if you are unsure about anything or are touching any files outside
# of your hk26-<team> directory.

# Checkout main and get your new changes
git checkout main
git pull
```

Feel free to ask if you need help with this.

## Getting help

If you run into problems with JASMIN accounts, you can access [JASMIN support](https://www.jasmin.ac.uk/help/contact/).
Make sure to put hk26 in the subject field. If you have problems with these instructions or the software, please email
[Mark Muetzelfeldt](mailto:mark.muetzelfeldt@reading.ac.uk), and put hk26 in the subject.
