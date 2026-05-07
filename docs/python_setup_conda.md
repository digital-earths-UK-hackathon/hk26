# Setting up a Python environment using `conda`

Log in to the notebook service: https://notebooks.jasmin.ac.uk/

Open a terminal:

Clone this repo:

```
git clone https://github.com/digital-earths-UK-hackathon/hk26.git
```

Create a Python env and add the kernel to the notebook service:

```
conda env create -f hk26_env.yaml
python -m ipykernel install --user --name hk26_env --display-name "Hackathon 2026 (conda)"
```

You can now create a new notebook with the "Hackathon 2026 (conda)" kernel to use all of the packages.

To install additional packages:

```
conda activate hk26_env
conda install <package_name>
```
