library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
nwk <- "newick_tree_modified.dist"
tree <- read.tree(nwk)