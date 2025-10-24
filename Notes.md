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