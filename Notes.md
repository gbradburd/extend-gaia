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




total accuracy & significance of accuracy gains (one histogram per sigma) 
total accuracy as function of sigma (mean & SE across sim reps per sigma)
accuracy as a function of time (one per sigma)
accuracy as a function of total extension (one per sigma)
accuracy of proportion correct extension (one per sigma)

log y axis in second one - to see more difference in error - still too close? make it into a table 
bold face the lowest one 

why decrease in error at deeper times? gaia thing but interesting to look at 

if wanna smash all graphs together in last one - make it one with all the lines color code and dash / dot them 




for the green
p values bigger in the one with those 
just write the labels on the left side or just on the bottom
p value bold and on the plot 
maybe legend function/argument?


        on indv subplot graphs fix the y axis so can see relationship across sigma 
        add log to the ones that logged
fix goofed graphg
        make sigma one a table 






add std dev on points - similar to adding another trace 
s       implify tabke more - 3 sig figs 
- boldface the one thats the best result 

could do a difference test - hypothesis testing that difference in error = 0 




sanity check accuracy of all the nodes - extended parent children etc etc 

mike does error in units of range size - explains why no relationhip there bc errpr re;ative to range size is smaller 

plots of total accuracy - mak sure nothing insane isnt happening why not relationship within sigma between depth of age and accuracy 

across simulations do absolute error rather than relative error - then see relationship b.w dispersal and error

look at accurcy through time and confirm that absolute error scales with diserpcal 
if does - then think we should thats good 
if dostn then issue somewhere in pipeline 


if absolute error fo all nodes doesnt increase with age - then skmethings wrong - worth mining finding whats wrong 

if id does - what do we do ? decide 





- run some absolute error on all nodes - see if it increases with age 
- 