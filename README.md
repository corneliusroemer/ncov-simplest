In order to run the analysis yourself, you need to `select all` on GISAID, `download` and then select `Input for the Augur pipeline`.

Put the downloaded `.tar` file into the `data` folder.

Then run Snakemake with `snakemake export -c4`.

The resulting `auspice.json` will be in the folder `auspice`.