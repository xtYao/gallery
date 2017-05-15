library(skitools)
library(gTrack)
library(gChain)
## library(gGnome)
load_all("~/gitLcl/gGnome/")

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
es = gg$edges
es[from==2 & to==4, cn:=1]
es[from==4 & to==2, cn:=2]
es[from==9 & to==7, cn:=1]
es[from==7 & to==9, cn:=2]
es[from==2 & to==3, cn:=3]
es[from==8 & to==7, cn:=3]
es[from==3 & to==4, cn:=3]
es[from==9 & to==8, cn:=3]
gg$hardset(segs=segs, es=es)

bg = bGraph$new(gg)
bg = bGraph$new(gg); gw = bg$walk() ## walk function is problematic!

## try new hurestic
G = bg$G
igraph::shortest_paths(G, 1, 5, weights = c(1e3, ))

## fake coverage
bins = gr.tile(gr, 200)
bins2 = bins %*% (segs %Q% (strand=="+"))
cov = sapply(bins2$cn, function(cn) rnorm(1, cn, sd=0.5))
bins2$cov = cov

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

## cov track
cov.td = gTrack(bins2, y.field = "cov")
cov.td$name = "cov ratio"
cov.td$circles = TRUE
cov.td$cex = 0.5
pdf("~/gitLcl/gallery/tDupAsDel.pdf", width=8, height=6, family="sans", fonts="Helvetica")
plot(c(cov.td, td), gr)
dev.off()

## now use gChain to plot the full history and change of haplotypes
help(gChain)
