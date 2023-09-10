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
ref_host_assignment = annot_df #targets::tar_read(tbl_refseq_annot)
rpm_table = targets::tar_read(rpm_table_virus_masked)
all_blastp = read_tsv("raw_data/all.blastp.txt")

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

tree_list = list()
plot_list = list()


for (clade in clades) {
  tree = read.tree(str_glue("raw_data/viral_genomes/phylogenetic_trees/{clade}.phy_phyml_tree.txt"))
  tree = midpoint.root(tree)
  tree_list[[clade]] = tree
  
  family_assignment = tibble(label = tree$tip.label) %>%
    left_join(select(all_blastp, qseqid, family_assignment), by=c("label"="qseqid")) %>%
    filter(!is.na(family_assignment)) %>%
    group_by(family_assignment) %>%
    group_modify(function(data, group) tibble(mrca=findMRCA_wrap(tree, data$label))) %>%
    ungroup()
  
  tree_with_meta = left_join(tree, left_join(tree_meta, 
                                             na.omit(select(all_blastp, label=qseqid, family_assignment)), 
                                             by="label"), 
                             by="label")
  
  g = ggtree(tree_with_meta, aes(color=invert_assoc)) + #family_assignment)) +
    geom_tiplab(aes(color=tag), 
                align=F, linetype = 1, linesize = 0.5, size=2.5, offset=0.1) +
    geom_tippoint(aes(shape=evidence_rpm), color="red", size=1.5, stroke = 0.8) +
    geom_treescale(width=0.5, offset=-10, fontsize=2, linesize=0.2)
  
  for (row in 1:nrow(family_assignment)) {
    n = family_assignment$mrca[row]
    f = family_assignment$family_assignment[row]
    g = g + geom_cladelab(node=n, label=f)
  }
  

  g = g +
    scale_shape_manual(values=c("rpm1k"=5, "rpm10k"=5, "none"=NULL)) +
    scale_color_manual(values=c(c("ref"="black",
                                  "core"="red",
                                  "uncertain"="purple",
                                  "non-core"="blue",
                                  "TRUE"="green",
                                  "FALSE"="#444444")),
                       na.value="black") +
    theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 0, unit = "mm"),
          legend.position = "none")
  
  g = g + xlim(c(0, layer_scales(g)$x$range$range[2] * 1.7))
  
  plot_list[[clade]] = g
  ggsave(filename=str_glue("output/SI_trees/{clade}.jpg"), plot=plot_list[[clade]],
         width = 180, height = 20 + 2 * Ntip(tree_list[[clade]]),
         units = "mm")
}

