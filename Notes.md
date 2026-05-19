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







look at statistics on distribution 
instead of doing difference - pro

do time results in R - if can get better time results in slow ones - if the slow as hell methods are more viable bc so extension 
rerun s1.1 0 or turn them from strigns into floats ...., bruh 
run more slimulations
figure otu where gaia does bad - seee how improves on those liek - bigger improvement?
keep running on more ones 
fix vpn issues 
figure out cluster more bc good to figure out
figure out python reticulate so can run tskit stuff in R
add extrea space to 0.2 file names so its consistent 




investigate examples with high error change between replicates - some examples dont change mucg bw replicates some do 
compare like proportions of nodes that were extended and by how much between examples both replicates and other sigmas 
threshold nodes that were extended the most and nodes with the highest difference in error 
    - did they have error in extension? how much extended? lots of parents not lots of parents? lots of children not lots of children? 
    how many trees did it get extended into?


- in s1.1 there is more variation between the args - some of them gaia does well on some it doesnt, when gaia does poorly extending helps ~ the same amount 
- trying to figure out why gaia is doing bad on these particular examples 
- to do with gaia not extension 

think more about how gaia actually works when were doing this 


something in the arg itself thats happening that we have to figure out 


get data on number of trees its in before and post extension 


try linear R
    - just doesnt work doesnt work 
    R: treeseq_sankoff_linear.c:341: plf_add: Assertion `fequals(slope[0], 0)' failed.
try in quadratic R
    - ets is slightly slower 
try in discrete R
    - 

try on fast gaia in argscape
    - ets is faster in fast gaia
try on quadratic gaia in argscape
    - ets is faster 
try on linear gaia in argscape 


try in python quadratic 
    - ets is slightly slower 
try in python linear
    - ets is slightly slower 
python discrete is not implemented 



in ets we add and remove less edges so why is gaia slower?

look into what sankoff era means ?
try linear on different arg 


- run t tests across simulations 
- why faster in argscape and not others?
    - dont really need to 
- mostly small reduction in speed 
- getting plots of total accuracy and significance of accuracy gains 
        - allow to zoom in on what the story is 
        - presumably accuracy does increase just more nuanced 
        - across simulations and sigma how does accuracy increase and for what time depths 
        - once know this be able to figure out sort of precicly what the story is 
- finish getting replications - go up from 9 to 20 replications 
- stick in python 