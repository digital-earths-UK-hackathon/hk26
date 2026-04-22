# Setting up a Python environment using `uv`

Log in to the notebook service: https://notebooks.jasmin.ac.uk/

Open a terminal:

Install `uv`

```
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Clone this repo:

```
git clone https://github.com/digital-earths-UK-hackathon/hk26.git
cd hk26
```

Create a Python env and add the kernel to the notebook service:

```
uv sync
uv run python -m ipykernel install --user --name hk26_env --display-name "Hackathon 2026 (uv)"
```

`uv sync` looks for `pyproject.toml` to see what Python pakcages need to be installed, and installs them.

You can now create a new notebook with the "Hackathon 2026 (uv)" kernel to use all of the packages.

If you need more packages, return to the `hk26` directory and run `uv add <package_name>`

If these steps do not work, you can try the [steps for installing using `conda`](python_setup_conda.md).
