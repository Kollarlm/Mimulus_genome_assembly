# This script is used to take the .delta file output as part of the nucmer software
# and make my own ggplot figures.

library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)
library(ggmagnify)
library(ggfx)


### Function for reading in .delta file 
### (adapted from https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-ggplot2/)

readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

## reading in .mdelta files (filtered .delta)
## Check to see if there is a way to adapt the function to read in anything with a common string
chr1 <- readDelta("./data/pseudochromosome_L1_S1_chr1.mdelta")
chr2 <- readDelta("./data/pseudochromosome_L1_S1_chr2.mdelta")
chr3 <- readDelta("./data/pseudochromosome_L1_S1_chr3.mdelta")
chr4 <- readDelta("./data/pseudochromosome_L1_S1_chr4.mdelta")
chr5 <- readDelta("./data/pseudochromosome_L1_S1_chr5.mdelta")
chr6 <- readDelta("./data/pseudochromosome_L1_S1_chr6.mdelta")
chr7 <- readDelta("./data/pseudochromosome_L1_S1_chr7.mdelta")
chr8 <- readDelta("./data/pseudochromosome_L1_S1_chr8.mdelta")
chr9 <- readDelta("./data/pseudochromosome_L1_S1_chr9.mdelta")
chr10 <- readDelta("./data/pseudochromosome_L1_S1_chr10.mdelta")
chr11 <- readDelta("./data/pseudochromosome_L1_S1_chr11.mdelta")
chr12 <- readDelta("./data/pseudochromosome_L1_S1_chr12.mdelta")
chr13 <- readDelta("./data/pseudochromosome_L1_S1_chr13.mdelta")
chr14 <- readDelta("./data/pseudochromosome_L1_S1_chr14.mdelta")

## I am skipping the filter step since I already filter with nucmer but here is the function 
## Adapted from the authors above.

# filterMum <- function(df, minl=1000, flanks=1e4){
#   coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
#     summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
#     ungroup %>% arrange(desc(rs)) %>%
#     mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
#   merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
#     mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
# }
# 
# mumgp.filt = filterMum(mumgp, minl=1e4)
# mumgp.filt %>% head %>% kable

## Plot

chr8_plot <- ggplot(chr8, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) + geom_segment() +
  geom_point(alpha=.5) + facet_grid(qid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=9),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_text(angle=0, size=9), axis.ticks.y=element_blank()) +
  xlab('LMC') + ylab('SWB') + scale_colour_brewer(palette='Set1')
ggmagnify(chr8_plot, ylim = c(4996241,7263731), xlim = c(926660,1449416),
          inset_xlim = c(10, 8000000), inset_ylim = c(10000000, 20000000), axes = FALSE,
          border = TRUE, shadow = TRUE)
