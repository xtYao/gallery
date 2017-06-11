library(skitools)
library(gUtils)
library(gTrack)
library(JaBbA)


## model chromothripsis
this.chrom = si2gr(hg_seqlengths())["1"]
gr = gaps(gr.sample(this.chrom, 20, len = 1)) %Q% (strand == '+') %&% this.chrom
#gr = gr[unique(c(1, which(width(gr)>1e6), length(gr)))] ##pick every other
gr = gr[unique(c(seq(1, length(gr), 2), length(gr)))] ##pick every other
gr = reduce(gr+2e6)

## shuffle strands
gr = gr[order(gr.stripstrand(gr))]
strand(gr) = ifelse(runif(length(gr))>0.5, '+', '-')
strand(gr)[c(1, length(gr))] = '+'
names(gr) = letters[1:length(gr)]

## random permutation starting with 1 and ending at length(gr)
walk = c(1, sample(2:(length(gr)-1)), length(gr))

## edges with formatting
edges = data.table(from = walk[-length(walk)], to = walk[-1])[,
        dist := abs(start(gr.start(gr[from], ignore.strand = FALSE))-start(gr.start(gr[to], ignore.strand = FALSE)))][
   ,":="(col = 'red', lwd = 2, h = 8, cex.arrow = 0, v = 1+0.1*log10(dist)*(0.2+runif(length(dist))),
         not.flat = TRUE, dangle.w = 20)]

adj = edges[, sparseMatrix(from, to, x = 1, dims = rep(length(gr), 2))]
p = all.paths(adj, all = TRUE)$paths[[1]]

title = paste0(ifelse(strand(gr)[p]=='+','', '-'), names(gr)[p], collapse = ' ')

## make gTrack of shattered chromosome
gt.chrom = gTrack(gr, y0 = 0, col = 'gray80', edges = edges, angle = 0, hadj.label = 0.5)

## compute coverage track of copy number data
cn = as(coverage(grbind(this.chrom, gr)), 'GRanges') %&% this.chrom
cov = gr.tile(cn, 1e4)
cov$score = rpois(length(cov), 40*cn$score[gr.match(cov, cn)])

gt.cov = gTrack(cov, 'score', y0 = 0)

ppdf(plot(c(gt.cov, gt.chrom), this.chrom), title = title)


#' Wednesday, Mar 15, 2017 07:02:56 PM

## model bfb cycles

## awkward code grab cytobands
bands = unlist(dat(karyogram())[[1]])
bands$arm = gsub('(.*[pq]).*', '\\1', bands$band)

## reduce to chromosome arms
arms = gr.reduce(bands, by = 'arm')
this.chrom = arms %Q% (arm == '1q') ## take 1q
gr = gaps(gr.sample(this.chrom, 6, len = 1)) %&% this.chrom ## randomly break

## color and make reciprocal ranges
strand(gr) = '+'
names(gr) = as.character(1:length(gr))
gr$col = brewer.master(length(gr), 'Oranges')
grf = gr.flipstrand(gr)
names(grf) = -(1:length(gr))
grf$col = brewer.master(length(gr), 'Greens')

chrom = c(gr, grf)
chrom$border = 'gray40'
chrom = gr.fix(chrom)

## hack to make a bfb string
my.bfb = split(chrom[c('1', '2', '3', '4', '5','-5', '-4', '-3',
              '3', '4', '-4', '-3', '3', '-3', '3', '4', '-4', '-3', '3', '4', '5', '-5',
              '-4', '-3', '-2')], '1')


## simulate poisson coverage
cn = as(coverage(grbind(arms[1:2], unlist(my.bfb))), 'GRanges')
cov = gr.tile(cn, 1e4)
cov$score = rpois(length(cov), 40*cn$score[gr.match(cov, cn)])
gt.cov = gTrack(cov, 'score', y0 = 0)
gt.cov$y1 = 500

## plot on reference genome
ppdf(plot(c(gt.cov, gTrack(my.bfb, draw.path = TRUE, lwd.border = 1, angle = 0)), '1'))


## make gChain to rearranged genome
spc = spChain(my.bfb)
my.bfb2 = links(spc)$y
values(my.bfb2) = values(my.bfb[[1]])[,c('col', 'border'), drop = FALSE]
names(my.bfb) = NULL
my.bfb2 = split(my.bfb2, '1')


## plot on rearranged genome
gt = gTrack(my.bfb2, draw.path = FALSE, lwd.border = 1,  angle = 0)
gt$draw.backbone = FALSE
ppdf(plot(gt, '1'), cex = c(0.5, 1))
