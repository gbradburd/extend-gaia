# Notes for Extend-Gaia!

## October

Current Useful Links:


### October 24th.

Created General Pipeline for generating Heat Maps off of ARGs constructed using GAIA and extend_haplotypes:

1. Run this to get tree sequence: [GAIA-Paper Simulations](https://github.com/blueraleigh/gaia-paper/tree/main/data/slim/continuous-space/uniform-landscape/gaussian-dispersal/simulations)
2. apply `.extend_haplotypes()` to get extended tree sequence (save both an extended tree-sequence and the original)
3. Run GAIA (python front end version from ARG-SCAPE) (download the tree-sequence .trees file and upload using tskit.load(path))
4. Generate Figures S3-S5 (Heat Kernels) [likely can find that code in GAIA-paper REPO] (do this for both the original tree sequence and the extended one so we can compare)
5. Celebrate because were going to get great results!!

## October 31st. Merry Spooky Season.

- Run Gaia (on data and extended) with default constants from the repo. Should we use default constants still or would there be a good reason to change them?

After some minimal difficulty we locally downloaded GAIA and got spatial data from the simulated ARGs by running GAIA on terminal from the instructions on the github repo. We decided to run spatial data using both `fast-gaia` and `gaia-quadratic` in case there are substantial differences in the methods and to compare. We also found how to run the heat kernel plots in the [manuscript repo](https://github.com/blueraleigh/gaia-paper/blob/223c7bd2449426ca580bc06ac7160c5b4a26ca19/data/slim/continuous-space/uniform-landscape/gaussian-dispersal/analysis/figs/ancestor-estimates.R). Plotting is in `R` so we will have to download R-studio!