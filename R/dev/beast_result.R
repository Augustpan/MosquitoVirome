library(tidyverse)
library(ggtree)
library(treeio)

options(ignore.negative.edge=TRUE)

lib_meta = read_csv("raw_data/lib_meta.csv")

aalb = read.beast("mosquito_phylogeography/aalb_mcc.txt")
asub = read.beast("mosquito_phylogeography/asub_mcc.txt")
asin = read.beast("mosquito_phylogeography/asin_mcc.txt")
ctri = read.beast("mosquito_phylogeography/ctri_mcc.txt")
cpip = read.beast("mosquito_phylogeography/cpip_mcc.txt")


asub = read.beast("mosquito_phylogeography/Aedes_albopictus_COX1.fa.aln.trees.mcc.txt")
labs = asub@phylo$tip.label
to_drop = labs[!labs %in% lib_meta$lib_id & labs != "Aedes_aegypti_COX1_NC_035159"]
asub_drop = drop.tip(asub, to_drop)
ggtree(asub_drop, mrsd = "2021-1-1") +
  geom_range("height_0.95_HPD", color='red', size=2, alpha=.3) + 
  theme_tree2() +
  labs(title="Armigeres subalbatus")

ggtree(asub) +
  geom_range("height_0.95_HPD", color='red', size=2, alpha=.3) + 
  theme_tree2() +
  labs(title="Armigeres subalbatus")