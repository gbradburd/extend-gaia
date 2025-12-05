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

## December 5th.
Gideon: 
Nomenclature:
 - T_og is the original tree that comes off the simulation, not simplified
 - T_og_s_e_true is the original tree that is simplified but without touching any of the nodes that show up in T_s (ie, it is the *true* extended simplified tree that T_s_e is an estimate of)
 - T_s is the simplified original tree, retaining only the modern-day sample nodes
 - T_s_e is the simplified original tree that is subsequently extended to contain unary nodes

Things to do and/or look at!

    - instead of side-by-side ancestry accuracy plots, look at a plot of difference in ancestor accuracy (only works if the nodes are identical btwn them)

    - look at accuracy plots for a specific subset of nodes:

            - just nodes that are extended in T_s_e

            - just nodes that are extended in T_s_e and those upstream and downstream by some delta

(maybe the two above but for T_og_s_e_true)

- How do we know which nodes are extended?  we go into T_s_e and identify all unary nodes. Any unary node has is one that has been extended!

    Can keep track of this by editing the function:

    def get_span_stats(ts, ets):
    time_map = {}
    added_span = np.zeros(ets.num_nodes)
    wrong_added_span = np.zeros(ets.num_nodes)
    for n in ts.nodes():
        # check times are unique except for samples
        assert n.time == 0.0 or n.time not in time_map
        time_map[n.time] = n.id
    for interval, t, et in ts.coiterate(ets):
        interval_length = interval[1] - interval[0]
        t_nodes = list(t.nodes())
        for n in et.nodes():
            if et.num_children(n) == 1:
                added_span[n] += interval_length
                # total_added_span += interval_length
            on = time_map[et.time(n)]
            if on not in t_nodes:
                assert et.num_children(n) == 1
                wrong_added_span[n] += interval_length
    return added_span, wrong_added_span
