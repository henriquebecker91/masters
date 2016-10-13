#!/usr/bin/R --vanilla
library(ggplot2)

n <- 100
wmax <- 1000
pmax <- 1000
s <- 42
set.seed(s)
w <- sample(1:wmax, n, replace=TRUE)
p <- sample(1:pmax, n, replace=TRUE)
t = data.frame(w = w, p = p)
pl = qplot(w, p, data = t, xlab = 'weight', ylab = 'profit', main = 'An uncorrelated instance')     
ggsave(filename = 'uncorrelated.eps', plot = pl, width = 210*0.75, height = 297*0.4, units = c('mm'), dpi = 300)

