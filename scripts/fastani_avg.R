#!/usr/bin/env Rscript

## Takes the raw output from FastANI and calculates average for each bidirectional pair of genomes
## Removes self-loops and multiple edges
## FastANI should have been run with the minFraction option included (0.85)

## $1: FastANI output file name
## $2: output file name (normal average ANI file)
## $3: output file name (ANI95 filtered, shifted for MCL)

args = commandArgs(trailingOnly=TRUE)

library(igraph)
library(dplyr)

data <- read.csv(args[1], sep="\t", header=F)

## Keep only first 3 columns and rename
data <- select(data, V1, V2, V3)
colnames(data) <- c("g1", "g2", "ANI")


## Replace .fasta in columns 1 and 2
data <- data.frame(lapply(data, function(x) { gsub("\\.fasta", "", x)}))

## Coerce ANI column to numeric
data$ANI <- as.numeric(as.character(data$ANI))

## Convert df into graph
G <- graph_from_data_frame(data, directed=F)


## check if graph is simple
## simple graphs are graphs which do not contain loops and multiple edges
## if not simple, then convert, getting mean of edges

#is_simple(G)
G.simple <- simplify(G, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb="mean")


## convert simple graph to df and write to file
newdata <- igraph::as_data_frame(G.simple, what="edges")
colnames(newdata) <- c("g1", "g2", "ANI")
write.table(newdata, file=args[2], sep="\t", quote=F, row.names=F)

## create file filtered to ANI>=95 for MCL clustering
newdata_shifted <- filter(newdata, ANI>=95) %>% transmute(g1, g2, ANI_shifted = (ANI-95)/100)
write.table(newdata_shifted, file=args[3], sep="\t", quote=F, row.names=F, col.names=F)
