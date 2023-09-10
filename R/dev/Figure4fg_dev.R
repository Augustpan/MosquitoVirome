library(tidyverse)
library(ComplexHeatmap)
library(circlize)

rpm_table = read_csv("output/rpm_table_virus_masked.csv")
vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))
lib_meta = read_csv("raw_data/lib_meta.csv")

cns_dist_pairwise = read_tsv("raw_data/viral_genomes_mapping/all_dist_p50.txt",
                             col_names=c("map_name1", "map_name2", "pident", "length")) %>%
  mutate(virus=str_extract(map_name1, "^(.+?)_.+?_len", group=1),
         lib1=str_extract(map_name1, "^.+?_(.+?)_len", group=1),
         lib2=str_extract(map_name2, "^.+?_(.+?)_len", group=1)) %>%
  select(-map_name1, -map_name2)


xy = SoDA::geoXY(lib_meta$latitude, lib_meta$longitude, unit=1000)
rownames(xy) = lib_meta$lib_id
geo_dist = dist(xy)
coi_tree = ape::read.tree("raw_data/non_contaminated_COI.trimal.aln.treefile")
coi_dist = cophenetic(coi_tree)
coi_dist = as.dist(coi_dist[lib_meta$lib_id,lib_meta$lib_id])
geo_dist_df = reshape2::melt(as.matrix(geo_dist), varnames = c("lib1", "lib2"), value.name = "geo_dist")
coi_dist_df = reshape2::melt(as.matrix(coi_dist), varnames = c("lib1", "lib2"), value.name = "coi_dist")

all_dist_df = cns_dist_pairwise %>%
  left_join(geo_dist_df) %>%
  left_join(coi_dist_df)

ds = sample_n(all_dist_df, 500000) %>% filter(length>2000)
g1 = ggplot(aes(x=geo_dist, y=pident), data=ds) + 
  ggrastr::rasterize(geom_point(size=1, alpha=0.1,  color="gray", shape=1), dpi=300) + 
  geom_smooth(method="lm", color="black", size=0.5) +  
#  scale_color_manual(values=c("TRUE"="black", "FALSE"="black")) +
  ylab("Percent genome identity\n(viruses)") +
  xlab("Spatial distance") +
  theme_bw()


coi_class = ifelse(ds$coi_dist>0.1, "s", "g")
g2 = ggplot(aes(x=coi_dist, y=pident, color=coi_class), data=ds) + 
  ggrastr::rasterize(geom_point(size=1, alpha=0.1, shape=1), dpi=300) + 
  geom_smooth(method="lm", size=0.5) + 
#  scale_color_manual(values=c("s"="black", "g"="black")) +
  ylab("Percent genome identity\n(viruses)") +

  xlab("Phylogenetic distance of mosquitos") +
  theme_bw() +
  theme(legend.position = "none")

f4f = g1
f4g = g2
g = cowplot::plot_grid(g2,g1, nrow=2)
ggsave(filename="output/Figure4fg.pdf", plot=g, width=3, heigh=4.6)
