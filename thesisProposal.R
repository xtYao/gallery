library(skitools)
library(skidb)
library(gTrack)
library(Flow)
load_all("~/git/gGnome")

## This is script for making the figures for thesis proposal
## Figure 1A. chromothripsis, chromoplexy, and BFB
## all chromothripsis loci in PCAWG
pcawgChromothripsis =
    dir("~/projects/PCAWG/Flow5/JaBbA.Nozzle/",
    pattern="chromothripsis.loci.rds",
    full.names=TRUE, recursive=TRUE)
