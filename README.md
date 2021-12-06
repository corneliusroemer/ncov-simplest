In order to run the analysis yourself, you need to `select all` on GISAID, `download` and then select `Input for the Augur pipeline`.

Put the downloaded `.tar` file into the `data` folder.

Set up your conda environment by `conda env create -f conda_env.yaml`

Then run Snakemake with `snakemake build -c0`.

The resulting `auspice.json` will be in the folder `auspice`.
