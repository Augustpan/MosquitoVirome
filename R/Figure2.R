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

draw_figure2 = function(virus_meta, 
                        ref_host_assignment,
                        rpm_table) {
  
  all_blastp = read_tsv("raw_data/all.blastp.txt")
  all_blastp$family_assignment = str_replace(all_blastp$family_assignment, "viridae$", "-")
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
    
    tree_with_meta = left_join(tree, tree_meta, by="label")
    
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
    
    g = ggtree(tree_with_meta, aes(color=invert_assoc), size=0.3) +
      geom_tiplab(aes(color=tag, label=""), 
                  align=T, linetype = 1, linesize = 0.3, size=1, offset=0.1) +
      geom_tippoint(aes(shape=evidence_rpm), color="red", size=0.6, stroke = 0.5) +
      geom_treescale(width=0.5, offset=-10, fontsize=2, linesize=0.2)
    for (row in 1:nrow(family_assignment)) {
      n = family_assignment$mrca[row]
      f = family_assignment$family_assignment[row]
      g = g + geom_cladelab(node=n, label=f, fontsize=1.5, align=T, lineheight=0.5)
    }
    g = g + scale_shape_manual(values=c("rpm1k"=5, "rpm10k"=5, "none"=NULL)) +
      scale_color_manual(values=c(c("ref"="#eeeeed", 
                                    "core"="red", 
                                    "uncertain"="purple",
                                    "non-core"="blue",
                                    "TRUE"="green",
                                    "FALSE"="#444444")),
                         na.value="black") +
      theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 0, unit = "mm"), 
            legend.position = "none")
    
    plot_list[[clade]] = g
  }
  
  #### merge plots ####
  
  library(cowplot)
  
  #### col1: picorna ####
  c1 = "picorna"
  
  # plots
  p1 = plot_list[[c1]]
  
  # trees
  t1 = tree_list[[c1]]
  
  p1_xmax = layer_scales(p1)$x$range$range[2]
  
  col1 = plot_grid(
    p1 + xlim(c(0, p1_xmax*1.5)) + theme(plot.margin = unit(c(0,0,0,0), "cm")), 
    labels = c("Picorna"),
    label_size = 6, vjust=0, hjust=0,
    ncol = 1,
    rel_heights = c(Ntip(t1))
  )
  
  #### col2: tombus & luteo ####
  c1 = "tombus"
  c2 = "luteo"
  
  # plots
  p1 = plot_list[[c1]]
  p2 = plot_list[[c2]]
  
  # trees
  t1 = tree_list[[c1]]
  t2 = tree_list[[c2]]
  
  p1_xmax = layer_scales(p1)$x$range$range[2]
  p2_xmax = layer_scales(p2)$x$range$range[2]
  
  col2 = plot_grid(
    p1 + xlim(c(0, p1_xmax*1.5)), 
    p2 + xlim(c(0, p2_xmax*1.5)), 
    labels = c("Tombus-Noda", "Luteo-Sobemo"),
    label_size = 6, vjust=0, hjust=0,
    ncol = 1,
    rel_heights = c(Ntip(t1), Ntip(t2))
  )
  
  #### col3: hepe & narna ####
  c1 = "hepe"
  c2 = "narna_v2"
  
  # plots
  p1 = plot_list[[c1]]
  p2 = plot_list[[c2]]
  
  # trees
  t1 = tree_list[[c1]]
  t2 = tree_list[[c2]]
  
  p1_xmax = layer_scales(p1)$x$range$range[2]
  p2_xmax = layer_scales(p2)$x$range$range[2]
  
  col3 = plot_grid(
    p1 + xlim(c(0, p1_xmax*1.5)), 
    p2 + xlim(c(0, p2_xmax*1.5)), 
    labels = c("Hepe-Virga", "Narna-Levi"),
    label_size = 6, vjust=0, hjust=0,
    ncol = 1,
    rel_heights = c(Ntip(t1), Ntip(t2))
  )
  
  #### col4: nido & flavi & astro & quenya & partiti ####
  c1 = "nido"
  c2 = "flavi"
  c3 = "astro"
  c5 = "partiti"
  
  # plots
  p1 = plot_list[[c1]]
  p2 = plot_list[[c2]]
  p3 = plot_list[[c3]]
  p5 = plot_list[[c5]]
  
  # trees
  t1 = tree_list[[c1]]
  t2 = tree_list[[c2]]
  t3 = tree_list[[c3]]
  t5 = tree_list[[c5]]
  
  p1_xmax = layer_scales(p1)$x$range$range[2]
  p2_xmax = layer_scales(p2)$x$range$range[2]
  p3_xmax = layer_scales(p3)$x$range$range[2]
  p5_xmax = layer_scales(p5)$x$range$range[2]
  
  col4 = plot_grid(
    p1 + xlim(c(0, p1_xmax*1.5)), 
    p2 + xlim(c(0, p2_xmax*1.5)), 
    p3 + xlim(c(0, p3_xmax*1.5)), 
    p5 + xlim(c(0, p5_xmax*1.5)), 
    labels = c("Nido", "Flavi", "Astro-Poty", "Partiti-Picobirna"),
    label_size = 6, vjust=0, hjust=0,
    ncol = 1,
    rel_heights = c(Ntip(t1), Ntip(t2), Ntip(t3), Ntip(t5))
  )
  
  #### col5: reo, birna & toti ####
  c1 = "reo"
  c2 = "birna"
  c3 = "toti_v2"
  
  
  # plots
  p1 = plot_list[[c1]]
  p2 = plot_list[[c2]]
  p3 = plot_list[[c3]]
  
  # trees
  t1 = tree_list[[c1]]
  t2 = tree_list[[c2]]
  t3 = tree_list[[c3]]
  
  p1_xmax = layer_scales(p1)$x$range$range[2]
  p2_xmax = layer_scales(p2)$x$range$range[2]
  p3_xmax = layer_scales(p3)$x$range$range[2]
  
  col5 = plot_grid(
    p1 + xlim(c(0, p1_xmax*1.5)), 
    p2 + xlim(c(0, p2_xmax*1.5)), 
    p3 + xlim(c(0, p3_xmax*1.5)), 
    labels = c("Reo", "Birna-Permutotetra", "Toti-Chryso"),
    label_size = 6, vjust=0, hjust=0,
    ncol = 1,
    rel_heights = c(Ntip(t1), Ntip(t2), Ntip(t3))
  )
  
  #### col6: bunya & orthomyxo ####
  c1 = "bunya"
  c2 = "orthomyxo"
  c3 = "parvo"
  # plots
  p1 = plot_list[[c1]]
  p2 = plot_list[[c2]]
  p3 = plot_list[[c3]]
  # trees
  t1 = tree_list[[c1]]
  t2 = tree_list[[c2]]
  t3 = tree_list[[c3]]
  p1_xmax = layer_scales(p1)$x$range$range[2]
  p2_xmax = layer_scales(p2)$x$range$range[2]
  p3_xmax = layer_scales(p3)$x$range$range[2]
  
  col6 = plot_grid(
    p1 + xlim(c(0, p1_xmax*1.5)), 
    p2 + xlim(c(0, p2_xmax*1.5)), 
    p3 + xlim(c(0, p3_xmax*1.5)), 
    labels = c("Bunya-Arena", "Orthomyxo", "Parvo"),
    label_size = 6, vjust=0, hjust=0,
    ncol = 1,
    rel_heights = c(Ntip(t1), Ntip(t2), Ntip(t3)-150)
  )
  
  #### col7: mono & yue ####
  c1 = "mono_v2"
  c2 = "yue"
  c3 = "baci-circo"
  c4 = "nudi"
  
  
  # plots
  p1 = plot_list[[c1]]
  p2 = plot_list[[c2]]
  p3 = plot_list[[c3]]
  p4 = plot_list[[c4]]
  
  # trees
  t1 = tree_list[[c1]]
  t2 = tree_list[[c2]]
  t3 = tree_list[[c3]]
  t4 = tree_list[[c4]]
  
  p1_xmax = layer_scales(p1)$x$range$range[2]
  p2_xmax = layer_scales(p2)$x$range$range[2]
  p3_xmax = layer_scales(p3)$x$range$range[2]
  p4_xmax = layer_scales(p4)$x$range$range[2]
  
  col7 = plot_grid(
    p1 + xlim(c(0, p1_xmax*1.5)), 
    p2 + xlim(c(0, p2_xmax*1.5)), 
    p3 + xlim(c(0, p3_xmax*1.5)), 
    p4 + xlim(c(0, p4_xmax*1.5)), 
    labels = c("Mono-Chu", "Qin-Yue", "Circo", "Nudi"),
    label_size = 6, vjust=0, hjust=0,
    ncol = 1,
    rel_heights = c(Ntip(t1), 40, Ntip(t3)-150, 40)
  )
  
  combined_plot = plot_grid(
    col1, col2, col3, col4, col5, col6, col7,
    ncol=7
  ) + theme(plot.margin = margin(t = 2, r = 2, b = 0, l = 2, unit = "mm"))
  
  combined_plot
}
