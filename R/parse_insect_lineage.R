library(tidyverse)
library(targets)
library(ggtree)
library(ape)
library(phytools)

findMRCA_wrap <- function(tree, tip) {
  if (length(tip) == 1)
    return(which(tree$tip.label == tip))
  else
    return(findMRCA(tree, tip))
}

select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1,
         tree$tip.label[element],
         tree$node.label[element-Ntip(tree)])
} 
select.tip.or.node = Vectorize(select.tip.or.node, vectorize.args = "element")

virus_meta = targets::tar_read(tbl_virus_meta)
ref_host_assignment = targets::tar_read(tbl_refseq_annot)
rpm_table = targets::tar_read(rpm_table_virus_masked)
lineage = read_csv("raw_data/insect_lineage.csv")

max_rpm = rpm_table %>%
  pivot_longer(cols=2:ncol(.), names_to="virus_name", values_to="rpm") %>%
  filter(rpm > 0) %>%
  group_by(virus_name) %>%
  summarise(max_rpm = max(rpm))
filter(max_rpm, max_rpm > 10000)
tree_meta = bind_rows(
  virus_meta %>%
    select(label = virus_name, is_core, uncertainty) %>%
    mutate(is_ref = F),
  ref_host_assignment %>%
    select(label, invert_assoc = is_core) %>%
    mutate(is_ref = T)) %>%
  replace_na(list(invert_assoc = F, is_core = F)) %>%
  mutate(rpm1k = label %in% filter(max_rpm, max_rpm > 1000)$virus_name,
         rpm10k = label %in% filter(max_rpm, max_rpm > 10000)$virus_name) %>%
  mutate(evidence_rpm = ifelse(rpm10k, 
                               "rpm10k",
                               ifelse(rpm1k, "rpm1k", "none"))) %>%
  mutate(tag = ifelse(is_ref, 
                      "ref", 
                      ifelse(is_core, 
                             "core", 
                             ifelse(is.na(uncertainty), 
                                    "non-core", 
                                    "uncertain"))))

clades = c("astro", "birna", "bunya", "flavi",
           "hepe", "luteo", "mono_v2", "narna_v2", 
           "nido", "orthomyxo", "partiti", "picorna",
           "reo", "tombus", "toti_v2", "yue", 
           "baci-circo", "nudi", "parvo")

annot_df = NULL
for (clade in clades) {
  tree = read.tree(str_glue("raw_data/viral_genomes/phylogenetic_trees/{clade}.phy_phyml_tree.txt"))
  tree = midpoint.root(tree)

  tree_with_meta = left_join(tree, tree_meta, by="label")
  cl = clade
  ldf = filter(lineage, clade==cl)
  if (nrow(ldf) > 1) {
    for (row in 1:nrow(ldf)) {
      mrca = findMRCA_wrap(tree, c(ldf$node1[row], ldf$node2[row]))
      tips = tree$tip.label[getDescendants(tree, mrca)]
      tips = tips[!is.na(tips)]
      annot_df = rbind(annot_df, tibble(label=tips, is_core=TRUE, evidence="MANUAL"))
    }
  }
}

annot_df = select(ref_host_assignment, label) %>%
  left_join(distinct(annot_df), by="label")
annot_df %>%
  write_csv("raw_data/refseq_host_annotation.csv")
