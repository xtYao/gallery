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

tDupNestedPath = c(1, 2, 3, 4, 2, 4, 2, 3, 4, 5)

gw = gWalks$new(segs = gg$segstats,
                paths = list(tDupNestedPath,
                             rev(map[as.character(tDupNestedPath)]),
                             1:5, 10:6),
                isCyc = rep(FALSE, 4),
                cn = rep(1, 4),
                str = rep(c("+", "-"), each=2))
## walk function is problematic! Problem in eclass!
td.gw = xtYao.gTrack(gTrack(gw$as.grl[c(1, 3)], draw.path=T, gr.labelfield = "tile.id", gr.colorfield = "tile.id"))


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
td$name = "nested\ntandem duplication"
td$yaxis = TRUE

## cov track
cov.td = gTrack(bins2, y.field = "cov")
cov.td$name = "cov ratio"
cov.td$circles = TRUE
cov.td$cex = 0.5
pdf("~/gitLcl/gallery/tDupAsDel.pdf", width=8, height=6, family="sans", fonts="Helvetica")
plot(c(td, td1), gr, legend.params = list(plot=FALSE))
dev.off()

## now use gChain to plot the full history and change of haplotypes
help(gChain)

##########################################
## a chromothripsis model
## based on 10 segments
bps.ct = sort(gr.sample(gr, 9, 1, replace=F))
## break into segments
gg1 = gGraph$new()$trim(gr)$addSegs(bps.ct)

## segments that remain
## remain = sample(1:10, 6)
remain = sample(seq(2, 10, 2), 5)
## orientation scramble
remainStr = sample(c("+", "-"), length(remain), replace=TRUE)
hB = gg1$hydrogenBonds()
map = hB[, c(setNames(from, to), setNames(to, from))]
## ## swap the scrambled segs
remainRc = map[as.character(remain)]
remain = ifelse(remainStr=="+", remain, remainRc)
remainRc = map[as.character(remain)]

nn = length(remain)
from = remain[1:(nn-1)]
to  = remain[2:nn]
segs1 = gg1$segstats
es1 = data.table(from=from, to=to, cn=1, type="aberrant", weight=width(segs1[from]))
## Attention: turn segment ends to breakpoints, notice the 1bp shift and flip strand
bp.from = gr.end(segs1[es1$from], ignore.strand=FALSE)[,c()]
bp.from[which(strand(bp.from)=="-")] = bp.from[which(strand(bp.from)=="-")] %-% 1
bp.from = gr.flipstrand(bp.from)
## the other end
bp.to = gr.start(segs1[es1$to], ignore.strand=FALSE)[,c()]
bp.to[which(strand(bp.to)=="+")] = bp.to[which(strand(bp.to)=="+")] %-% 1
juncs1 = grl.pivot(GRangesList(bp.from, bp.to))

gg1$addJuncs(juncs1)
segs1$cn = ifelse(seq_along(segs1) %in% c(remain, remainRc), 2, 1)
segs1

es1 = gg1$edges
es1[type=="reference", cn := 1]
es1

gg1$hardset(segs=segs1, es=es1)
saveRDS(gg1, "chromothripsis.gg.rds")
plot(gg1) ## good to go
gg1$layout()

hB1 = gg1$hydrogenBonds()
map1 = hB1[, c(setNames(from, to), setNames(to, from))]

gw1 = gWalks$new(segs = gg1$segstats,
                 paths = list(remain,
                              rev(map1[as.character(remain)]),
                              1:10, 20:11),
                 isCyc = rep(FALSE, 4),
                 cn = rep(1, 4),
                 str = rep(c("+", "-"), each=2))

td1.gw = xtYao.gTrack(gTrack(gw1$as.grl[c(1, 3)], draw.path=T,gr.colorfield="tile.id", gr.labelfield="tile.id"))

td1 = gg1$td

td1$y0 = 0
td1@data[[1]]$segment  = rep(LETTERS[1:(length(segs1)/2)], 2)
td1$gr.colorfield = "segment"
td1@data[[1]]$lbl = rep(LETTERS[1:(length(segs1)/2)], 2)
td1$xaxis.ticklen = 0.25
td1$xaxis.chronly = T
td1$xaxis.unit = 1e6
td1$xaxis.suffix = "Mb"
td1$xaxis.round = 2
td1$xaxis.interval = 5e5
td1$xaxis.cex.tick = 0.75
td1$xaxis.cex.label = 0
td1$yaxis.pretty = 2
td1$yaxis.cex = 0.75
td1$lwd.border = 2
td1$gr.labelfield = "lbl"
td1$gr.cex.label = 1.2
td1$sep.lwd = 0.5
td1$name = "chromothripsis"

## TODO: the scale of event!!!
plot(td1, gr)

#####################################
## chromoplexy
bp2 = sort(gr.sample(gr, 3, 1, FALSE))
## junctions
bps2 = bp2[c(1,2,2,3,3,1)]
strand(bps2) = c("-", "+", "-", "+", "-", "+")
juncs2 = split(bps2, rep(1:3, each=2))
## make graph
gg2 = gGraph$new()$trim(gr)$addJuncs(juncs2)
plot(gg2)

chromoplexy.ix = c(1, 3, 2, 4)
hB2 = gg2$hydrogenBonds()
map2 = hB2[, c(setNames(from, to), setNames(to, from))]

gw2 = gWalks$new(segs = gg2$segstats,
                 paths = list(chromoplexy.ix,
                              rev(map2[as.character(chromoplexy.ix)]),
                              1:4, 8:5),
                 isCyc = rep(FALSE, 4),
                 cn = rep(1, 4),
                 str = rep(c("+", "-"), each=2))
td2.gw = xtYao.gTrack(gTrack(gw2$as.grl[c(1, 3)], draw.path=T, gr.colorfield="tile.id", gr.labelfield="tile.id"))

## make figure
td2 = gg2$td

td2$y0 = 0
td2@data[[1]]$segment  = rep(LETTERS[1:4], 2)
td2$gr.colorfield = "segment"
td2@data[[1]]$lbl = rep(LETTERS[1:4], 2)
td2$xaxis.ticklen = 0.25
td2$xaxis.chronly = T
td2$xaxis.unit = 1e6
td2$xaxis.suffix = "Mb"
td2$xaxis.round = 2
td2$xaxis.interval = 5e5
td2$xaxis.cex.tick = 0.75
td2$xaxis.cex.label = 0
td2$yaxis.pretty = 2
td2$yaxis.cex = 0.75
td2$lwd.border = 2
td2$gr.labelfield = "lbl"
td2$gr.cex.label = 1.2
td2$sep.lwd = 0.5
td2$name = "chromoplexy"

####################################
## BFB
## may have to rewrite gWalks class
bfb.ix = as.numeric(c('1', '2', '3', '4', '5','-5', '-4', '-3',
           '3', '4', '-4', '-3', '3', '-3', '3', '4',
           '-4', '-3', '3', '4', '5', '-5','-4', '-3', '-2'))

segs3 = gg$segstats
hB3 = gg$hydrogenBonds()
map3 = hB3[, c(setNames(from, to), setNames(to, from))]
bfb.ix = ifelse(bfb.ix>0, bfb.ix, map3[as.character(abs(bfb.ix))])

bfb.gw = gWalks$new(segs = segs3,
                    paths = list(bfb.ix, rev(map[as.character(bfb.ix)]), 1:5, 10:6),
                    isCyc = rep(FALSE, 4),
                    cn = rep(1,4))

gg3 = bfb.gw$gw2gg()
bfb.grl = bfb.gw$gw2grl()

td3.gw = xtYao.gTrack(gTrack(bfb.grl[c(1, 3)], draw.path=T, gr.labelfield = "tile.id", gr.colorfield="tile.id"))

td3 = gg3$td

td3$y0 = 0
td3@data[[1]]$segment  = rep(LETTERS[1:5], 2)
td3$gr.colorfield = "segment"
td3@data[[1]]$lbl = rep(LETTERS[1:5], 2)
td3$xaxis.ticklen = 0.25
td3$xaxis.chronly = T
td3$xaxis.unit = 1e6
td3$xaxis.suffix = "Mb"
td3$xaxis.round = 2
td3$xaxis.interval = 5e5
td3$xaxis.cex.tick = 0.75
td3$xaxis.cex.label = 0
td3$yaxis.pretty = 2
td3$yaxis.cex = 0.75
td3$lwd.border = 2
td3$gr.labelfield = "lbl"
td3$gr.cex.label = 1.2
td3$sep.lwd = 0.5
td3$name = "breakage\nfusion\nbridge"

pdf("~/gitLcl/gallery/tDupNested.pdf",
    width=12, height=12, family="sans", fonts="Helvetica")
plot(c(td, td.gw), gr, legend.params = list(plot=FALSE))
dev.off()

pdf("~/gitLcl/gallery/chromothripsis.pdf",
    width=12, height=12, family="sans", fonts="Helvetica")
plot(c(td1, td1.gw), gr, name=NULL, legend.params = list(plot=FALSE))
dev.off()

pdf("~/gitLcl/gallery/chromoplexy.pdf",
    width=12, height=12, family="sans", fonts="Helvetica")
plot(c(td2, td2.gw), gr, legend.params = list(plot=FALSE))
dev.off()

pdf("~/gitLcl/gallery/bfb.pdf",
    width=12, height=12, family="sans", fonts="Helvetica")
td3$name=NULL
plot(c(td3, td3.gw), gr, legend.params = list(plot=FALSE))
dev.off()

pdf("~/gitLcl/gallery/complex.pdf",
    width=12, height=18, family="sans", fonts="Helvetica")
plot(c(td, td.gw, td1, td1.gw, td2, td2.gw, td3, td3.gw), gr, legend.params = list(plot=FALSE))
dev.off()

xtYao.gTrack = function(td){
    td$xaxis.ticklen = 0.25
    td$xaxis.chronly = T
    td$xaxis.unit = 1e6
    td$xaxis.suffix = "Mb"
    td$xaxis.round = 2
    td$xaxis.interval = 1e6
    td$xaxis.cex.tick = 0.75
    td$xaxis.cex.label = 0
    td$yaxis.pretty = 2
    td$yaxis.cex = 0.75
    td$lwd.border = 2
    td$gr.cex.label = 1.2
    td$sep.lwd = 0.5
    return(td)
}

## ## remove the copy number of removed segments
## segs1$cn[c(remain, remainRc)] = 2
## ## remove the copy number of old ref edges
## es1 = gg1$edges
## es1$cn=1
## nn = length(remain)

## from = c(remain[1:(nn-1)], remainRc[(nn-1):1])
## to = c(remain[2:nn], remainRc[nn:2])
## cn = 1
## type = "aberrant"
## weight = width(segs1)[from]
## es1 = rbind(es1, data.table(from=from, to=to, cn=cn, type=type, weight=weight))
## G1 = igraph::make_graph(t(as.matrix(es1[, .(from, to)])))
## bp1 = es1[type=="aberrant", gr.end(segs1[from], ignore.strand=F)]
## ## wait!!! need to shift - nodes 1bp left
## ## bp1 = ifelse(strand(bp1)=="-", bp1 %-% 1, bp1) ## how come this is breaking the mem??
## bp1[which(strand(bp1) == "-")] = bp1[which(strand(bp1) == "-")] %-% 1
## bp2 = es1[type=="aberrant", gr.start(segs1[to], ignore.strand=F)]
## bp2[which(strand(bp2) == "-")] = bp2[which(strand(bp2) == "-")] %-% 1
## junc1 = grl.pivot(GRangesList(list(bp1[,c()], bp2[,c()]))) ## TODO: this is duplicated!!!

## let's do it from junctions







