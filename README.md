# Repository for the UK km-scale hackathon 2026

## Getting started on JASMIN

If you haven't already set up Python on JASMIN, see [this guide](docs/python_setup_conda.md).

If planning to commit any code to community repositories on GitHub, you may want to check [this GitHub setup guide](docs/github_setup.md)

## Dataset catalog

You can then open the datasets contained in [this intake 0.7 catalog](https://digital-earths-global-hackathon.github.io/catalog/).
You can find all Unified Model by looking for entried that start with `um_`, and those that were created for this hackathon
by their `_hk26` suffix.

## K-Scale R&D datasets availability and access

For the DYAMOND-3 year long simulations the native .pp outputs from the Met Office Unified Model are currently available on JASMIN via the kscale group workspace. If you are not already part of the group

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
