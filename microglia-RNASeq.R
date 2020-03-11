library(edgeR)
library(rgr)

file <- "/Users/bicakm01/Desktop/microgliaData/all.final"

#counts <- read.delim(file)

data <- read.delim(file)
rowNames <- as.character((data$GeneName))
data <- data.matrix(data)
row.names(data) <- rowNames
data <- data[,-1]
dim(data)

tmp.data <- remove.na(data)
x <- tmp.data$x[1:tmp.data$n]

data2 <- remove.na(data)
xx <- c(15, 39, 18, 16, NA, 53)
temp.x <- remove.na(xx)
x <- temp.x$x[1:temp.x$n]

## to recover the other values returned
n <- temp.x$n
m <- temp.x$m

#the number of rows removed from xx.
nna <- temp.x$nna


d0 <- DGEList(data2)
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)

snames <- colnames(counts)

samples <- substr(snames, 1, nchar(snames) -2)
time <- substr(snames, nchar(snames) -1, nchar(snames) -1)

samples
time

group <- interaction(samples, time)

group

plotMDS(d, col = as.numeric(group))

mm <- model.matrix