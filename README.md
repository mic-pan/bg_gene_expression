# Supplementary code for "Thermodynamically-consistent, reduced models of gene regulatory networks"

## Scripts
This repository contains the Julia code required to run the bond graph model in the study *Thermodynamically-consistent, reduced models of gene regulatory networks*. The following files will generate the results figures as indicated below:
| Script                           | Figures         |
|----------------------------------|-----------------|
|`elongation_model_comparison.jl`  | Fig. 7          |
|`translation_model_comparison.jl` | Figs. 8, 9, S2  |
|`gene_regulation_sims.jl`         | Figs. 10-11, S3 |
|`resource_dependence.jl`          | Figs. 12-14     |
|`heterogeneity.jl`                | Fig. 15         |
|`toggle_phase_plane.jl`           | Fig. S1         |


## Instructions for running code
The code can be run using a Julia installation, using the steps below. The code has been tested on Julia 1.9.1.

1. Install Julia (https://julialang.org/downloads/) and Git (https://git-scm.com/downloads)
2. Clone this repository locally by running the following in the terminal: `git clone https://github.com/mic-pan/bg_gene_expression`
3. Start Julia: `julia`
4. Change into the code directory: `cd("bg_gene_expression")`
5. Press the `]` key to enter the package manager
6. Activate the environment: `activate env`
7. Install the required packages: `instantiate`
8. Exit the package manager by pressing the backspace key
9. The scripts can be run by using the command `include("script_name.jl")`
