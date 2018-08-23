# o2r examples

Examples for [Executable Research Compendia](https://o2r.info/erc-spec) and compatible workspaces.
Examples are ZIP archives to allow direct upload to the [o2r reproducibility service](https://o2r.info/architecture/), or directories which must be zipped before upload.

The o2r demo server is available at http://o2r.uni-muenster.de/.

For instructions on how to run the o2r reference implementation locally, see https://github.com/o2r-project/reference-implementation.

**NOTE:** Licenses of examples vary, please check carefully if they suit your needs.

We want to thank everyone who contributed a workspace! Why not [add your own?](https://o2r.info/almost/)

## Workspaces

### Journal articles or short papers

- [Garcia, M., Portney, K., and Islam, S.: A question driven socio-hydrological modeling process, Hydrol. Earth Syst. Sci., 20, 73-92, https://doi.org/10.5194/hess-20-73-2016, 2016.](workspaces/Aquestiondrivenprocess) 
[License: [Creative Commons Attribution 3.0 License](https://creativecommons.org/licenses/by/3.0/)]
- [Fraza, E., Elsner, J. B., and Jagger, T. H.: A space–time statistical climate model for hurricane intensification in the North Atlantic basin, Adv. Stat. Clim. Meteorol. Oceanogr., 2, 105-114, https://doi.org/10.5194/ascmo-2-105-2016, 2016.](workspaces/Aspacetimemodel) [License: [Creative Commons Attribution 3.0 License](https://creativecommons.org/licenses/by/3.0/)]
- [Nüst, Daniel. (2018, January 8). Open environmental data analysis with senseBox, openSenseMap, Jupyter Notebook, RStudio, and BinderHub. Zenodo.](workspaces/sensebox-opensensemap) is included in a shorted version, see [Zenodo](http://doi.org/10.5281/zenodo.1139929) [License: [Creative Commons Attribution 4.0 License](https://creativecommons.org/licenses/by/4.0/)]

### Demos

#### R Markdown

- [Minimal workspace as ZIP archive](workspaces/minimal-rmd.zip) with two files: `main.Rmd` in R Markdown format as the main analysis, and `display.html` in HTML for display [License: [Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)].
- [R Markdown workspace with data file](workspaces/rmd-data) with three files: `main.Rmd` in R Markdown format as the main analysis, `data.csv` providing the raw data for the analysis, and `display.html` in HTML for display [License: [Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)]. The workspace is also available in three variants to demonstrate the capabilities of the platform:
  - [R Markdown workspace with randomised plot](workspaces/rmd-data-random) changes the order of the included bar plot. This example demonstrates differences between the uploaded display file and the one created by the platform in a wrongly designed workflow (i.e. using random effects in a way breaking reproducibility).
  - [R Markdown workspace with wrong display file](workspaces/rmd-data_wrong-displayfile) has a mismatch between the plot created in the code and the plot included in the display file. This example demonstrates differences between the uploaded display file and the one created by the platform.
  - [R Markdown workspace with different dataset](workspaces/rmd-data-other) is a direct copy of the original R Markdown workspace but with different data values in `data.csv`. This example can be used to demonstrate a simple substitution of tabular data from a plain text format.
- [R Markdown workspace with `RData` files](workspaces/rmd-rdata) is a minimal workspace to demonstrate [o2r-inspecter](https://github.com/o2r-project/o2r-inspecter), exposing contents of the `RData` format in the o2r platform UI.
- [R Markdown workspace with text differences](workspaces/rmd-textdiff) is a minimal workspace to demonstrate the text difference checking.

#### R

- [Minimal workspace as ZIP archive](workspaces/minimal-script.zip) with two files, `main.R` and `display.png`, with data from http://www.budgetshippingcontainers.co.uk [License: [Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)]

#### Networking

The ERC's runtime container is a standalone environment, and it must not rely on any online resources to work.
The following two are a minimal demonstrator for this conceptual decision.

- [ping workspace](workspaces/ping) contains a `Dockerfile` (so the reproducibility service will not create one) which runs the [`ping`](https://en.wikipedia.org/wiki/Ping_(networking_utility)) tool to [`127.0.0.1`](https://en.wikipedia.org/wiki/Localhost), i.e. the localhost loopback, for 60 seconds. The other included files are merely there to make the metadata extraction work. This workspace can be successfully executed with the reproducibility service. [License: [Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)]
- [ping workspace "bad"](workspaces/ping-bad) contains files just as the previous workspace, but tries to ping an online URL (`o2r.info`). This workspace results in a failed execution with the reproducibility service. [License: [Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)]

## ERC

- [metadata test](ERC/metadata.zip) (ready to use ZIP archive) [License: [Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)]
- [ping test](ERC/ping.zip) (ready to use ZIP archive) [License: [Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)]

## Tools

The directory `corpus` contains an R Markdown document to upload the full corpus of demo papers to a local or remote instance of the o2r reference implementation.
