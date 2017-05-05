library(skitools)
library(gTrack)
library(gChain)
library(gGnome)

## TODO: write this creation of junctions into a function
## get a 1Mb gr
gr = gr.rand(1e6, GENOME)

## get 4 breakpoints in gr
bps = gr.sample(gr, 4, wid=1, replace=F)
bpsPlus = bps
strand(bpsPlus) = "+"
bpsMinus = gr.flipstrand(bpsPlus)

## get segmentation
segs = gr.breaks(gr, bps)
segs$cn=c(2, 4, 3, 4, 2)
strand(segs) = "+"
segs = c(segs, gr.flipstrand(segs))

## get junctions
## 1+,4-; 2-,3+
juncs = GRangesList(
    c(bpsPlus[1], bpsMinus[4]),
    c(bpsMinus[2], bpsPlus[3])
)

## create gGraph
