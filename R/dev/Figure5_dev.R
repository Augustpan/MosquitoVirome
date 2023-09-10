library(tidyverse)
library(ape)
library(ggtree)

lib_meta = read_csv("output/lib_meta_merged.csv")

rpm_table = read_csv("output/rpm_table_virus_masked.csv")
vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))

cst_profile = read_csv("output/cst_all.csv")

xy = SoDA::geoXY(lib_meta$latitude, lib_meta$longitude, unit=1000)
rownames(xy) = lib_meta$lib_id
spat_dist = as.matrix(dist(xy))

#cst_profile = read_csv("output/virus_meta_merged.csv")
cst_profile = filter(vm, is_core)
color_mapping = setNames(
  c("#3F5388", "#4DBBD6", "#E64B37", "#019F87", "#F29B7F", "#F29B7F", "#F29B7F", "#F29B7F", "#F29B7F"),
  c("Armigeres", "Anopheles", "Culex", "Aedes", "Mansonia", "Mimomyia", "Coquillettidia", "Lutzia", "Others")
)

if (F) {
  for (v in unique(cst_profile$virus_name)) {
    tmp = rpm_table %>%
      select(lib_id, v=v) %>%
      filter(v>0) %>%
      left_join(select(lib_meta, lib_id, genus, species)) %>%
      mutate(genus = ifelse(genus %in% c("Mansonia", "Mimomyia", "Coquillettidia", "Lutzia"), "Others", genus)) %>%
      arrange(genus, species, desc(v)) %>%
      mutate(rank=1:nrow(.))
    
    g = ggplot(aes(x=rank, y=log10(v+1), fill=genus), data=tmp) + 
      geom_bar(stat="identity") +
      xlab("Individual mosquito") +
      ylab("log(RPM+1)") +
      scale_fill_manual(values=color_mapping) + 
      theme(legend.position = "none")
    ggsave(filename=str_glue("output/virus_strain_abundance/cst_{v}.pdf"), plot=g, width=3, height=3)
  }
}

mantel_df = NULL
results = NULL
mos_coi_tree = ape::read.tree("raw_data/non_contaminated_COI.trimal.aln.treefile")

for (v in unique(cst_profile$virus_name)) {

  tree_fname = str_glue("raw_data/viral_genomes_mapping/split_by_virus_p80/{v}.phy_phyml_tree.txt")
  if (!file.exists(tree_fname)) {
    print(v)
    next
  }
  try({
    tree = tree_fname %>%
      ape::read.tree() %>%
      phytools::midpoint.root()
    
    #sp = str_split(tree$tip.label, "_", simplify = T)
    #tree = ape::drop.tip(tree, tree$tip.label[which(sp[,1]=="REF")])
    sp = str_split(tree$tip.label, "_", simplify = T)
    tree$tip.label = sp[,2]
    htree = drop.tip(mos_coi_tree, setdiff(mos_coi_tree$tip.label, tree$tip.label))
    df = rpm_table %>%
      select(lib_id, v=all_of(v)) %>%
      # filter(v>0) %>%
      left_join(select(lib_meta, lib_id, genus, species)) %>%
      mutate(genus = ifelse(genus %in% c("Mansonia", "Mimomyia", "Coquillettidia", "Lutzia"), "Others", genus)) %>%
      arrange(genus, species, desc(v)) %>%
      mutate(rank=1:nrow(.)) %>%
      mutate(label=lib_id)
    
    mos = top5[5]
    to_drop = setdiff(tree$tip.label, intersect(filter(lib_meta, species==mos)$lib_id, tree$tip.label))
    tree = drop.tip(tree, to_drop)
#    try({
      vh_links = table(tree$tip.label, tree$tip.label)
      vdist = cophenetic(tree)
      hdist = cophenetic(htree)
      sdist = as.dist(spat_dist[tree$tip.label,tree$tip.label])
      nm = rownames(vdist)
      vdist = as.dist(vdist[nm,nm])
      hdist = as.dist(hdist[nm,nm])
      gvh = ggplot() + 
        geom_point(aes(x=as.vector(hdist), y=as.vector(vdist)), size=0.75, alpha=0.3) +
        geom_smooth(aes(x=as.vector(hdist), y=as.vector(vdist)), method="lm", size=.75) + 
        xlab("Host distance") +
        ylab("Virus distance") +
        theme(axis.title = element_text(size=6),
              axis.text = element_text(size=5))
      gvs = ggplot() + 
        geom_point(aes(x=as.vector(sdist), y=as.vector(vdist)), size=0.75, alpha=0.3) +
        geom_smooth(aes(x=as.vector(sdist), y=as.vector(vdist)), method="lm", size=.75) +
        xlab("Spatial distance") +
        ylab("Virus distance") +
        theme(axis.title = element_text(size=6),
              axis.text = element_text(size=5))
      ghs = ggplot() + 
        geom_point(aes(x=as.vector(sdist), y=as.vector(hdist)), size=0.75, alpha=0.3) +
        geom_smooth(aes(x=as.vector(sdist), y=as.vector(hdist)), method="lm", size=.75) +
        xlab("Spatial distance") +
        ylab("Host distance") +
        theme(axis.title = element_text(size=6),
              axis.text = element_text(size=5))
      
      gm = cowplot::plot_grid(gvh, gvs, ghs, nrow=1)
      ggsave(filename = str_glue("output/virus_dist_plots/{v}.jpg"),
             plot = gm,
             width = 120,
             height = 40,
             units = "mm")
      x = vegan::mantel.partial(vdist, hdist, sdist)
      xs = vegan::mantel.partial(vdist, sdist, hdist)
      xh = vegan::mantel(hdist, sdist)
      mantel_df = rbind(mantel_df, c(v, 
                                     x$statistic, x$signif, 
                                     xs$statistic, xs$signif,
                                     xh$statistic, xh$signif
                                     ))
      #D <- prepare_paco_data(hdist, vdist, vh_links)
      #D <- add_pcoord(D)
      #D <- PACo(D, nperm=100, seed=42, method="r0")
      
      #results = rbind(results, c(v, D$gof$p, D$gof$ss, D$gof$n))
#    })
    
    df = left_join(tree, df, by="label")
    g = ggtree::ggtree(df) + 
      ggtree::geom_tippoint(aes(color=genus), size=3) +
      ggtree::geom_treescale() +
      scale_color_manual(values=color_mapping)
    ggsave(filename=str_glue("output/virus_strain_tree_genus/cst_tree_{v}.pdf"), plot=g, width=3, height=3)
  })
}
med_rpm = rpm_table %>% 
  pivot_longer(cols=2:ncol(.), names_to="virus_name") %>% 
  filter(value>0) %>% 
  group_by(virus_name) %>% 
  summarise(med_lg_rpm = median(log10(value+1)))

mantel_df = mantel_df %>% as.data.frame() %>% as_tibble()
mantel_df$V2 = as.double(mantel_df$V2)
mantel_df$V3 = as.double(mantel_df$V3)
mantel_df$V4 = as.double(mantel_df$V4)
mantel_df$V5 = as.double(mantel_df$V5)
mantel_df$V6 = as.double(mantel_df$V6)
mantel_df$V7 = as.double(mantel_df$V7)
mantel_df = rename(mantel_df, virus_name=V1, r=V2, p=V3, rs=V4, ps=V5, rh=V6, ph=V7)

tmp =left_join(mantel_df, med_rpm)
g=ggplot() + 
  geom_point(aes(x=r,y=log10(p), color=r),size=3,data=mantel_df) + 
  geom_hline(yintercept = log10(0.05), linetype=2) +
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab") + 
  theme_bw() +
  xlab("Mantel Correlation") +
  ylab("log(P-value)")
ggsave(filename="output/vhcorr.pdf", plot=g, height=3, width=4.2)

g=ggplot() + 
  geom_point(aes(x=rs,y=log10(ps), color=rs),size=3,data=mantel_df) + 
  geom_hline(yintercept = log10(0.05), linetype=2) +
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab") + 
  theme_bw() +
  xlab("Mantel Correlation") +
  ylab("log(P-value)")
ggsave(filename="output/vscorr.pdf", plot=g, height=3, width=4.2)

g=ggplot(aes(x=r, y=med_lg_rpm), data=tmp) +
  geom_point(aes(color=r), size=3) + 
  geom_smooth(method="lm") +
  xlab("Mantel Corrlation")+
  ylab("Median log(RPM+1)") + 
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab") +
  theme_bw()
ggsave(filename="output/rpmcorr.pdf", plot=g, height=3, width=4.2)

mantel_df %>% write_csv("output/mantel_corr.csv")