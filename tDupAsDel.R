library(skitools)
library(gTrack)
library(gChain)
library(gGnome)

## TODO: write this creation of junctions into a function
## get a 1Mb gr
set.seed(123)
gr = gr.rand(1e6, GENOME)
## get 4 breakpoints in gr
bps = sort(gr.sample(gr, 4, 1, replace=F))
## get junctions
## 1+,4-; 2-,3+
bpix = c(1, 4,
         2, 3)
bpor = c("+", "-",
         "-", "+")
juncs = bps[bpix]
strand(juncs) = bpor
juncs = split(juncs, f=rep(1:2, each=2))
## create gGraph
gg = gGraph$new()$trim(gr)$addJuncs(juncs)
## modify copy numbers
segs = gg$segstats
segs$cn = rep(c(2,4,3,4,2), 2)

## fine tune the plot
td = gg$td
td@data[[1]]$segment  = rep(LETTERS[1:5], 2)
td$gr.colorfield = "segment"
td@data[[1]]$lbl = rep(LETTERS[1:5], 2)
td$xaxis.ticklen = 0.25
td$xaxis.chronly = T
td$xaxis.unit = 1e6
td$xaxis.suffix = "Mb"
td$xaxis.round = 2
td$xaxis.interval = 5e5
td$xaxis.cex.tick = 0.75
td$yaxis.pretty = 2
td$yaxis.cex = 0.75
td$lwd.border = 2
td$gr.labelfield = "lbl"
td$gr.cex.label = 1.2
td$sep.lwd = 0.5
pdf("~/gitLcl/gallery/tDupAsDel.pdf", width=8, height=6, family="sans", fonts="Helvetica")
plot(td, gr)
dev.off()
