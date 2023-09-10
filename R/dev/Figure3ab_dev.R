library(tidyverse)

# load raw data
arbovirus = read_csv("raw_data/arbovirus.csv")
lib_meta = read_csv("raw_data/lib_meta.csv")
rpm_table = read_csv("output/rpm_table_virus_masked.csv")
read_count = read_tsv("raw_data/norRNA_reads.txt")
clim_data = read_csv("raw_data/climate.csv")

vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))

# reordering datasets
rpm_table = arrange(rpm_table, lib_id)
all_meta = lib_meta %>%
  left_join(read_count) %>%
  left_join(clim_data) %>%
  arrange(lib_id)

# should be exactly in the same order
assertthat::are_equal(rpm_table$lib_id, all_meta$lib_id)

# cluster individual samples within 100km as one site
#   hclust: average linkage (UPGMA)
#   cutree: 100km
geo_cluster = SoDA::geoXY(all_meta$latitude, all_meta$longitude, unit = 1000) %>%
  dist() %>%
  hclust(method = "average") %>%
  cutree(h=100)

# assign cluster info to each sample
all_meta$geo_cluster = paste0("c", geo_cluster)

arbo_table = select(rpm_table, lib_id, any_of(arbovirus$virus_name)) %>% 
  pivot_longer(cols=2:ncol(.)) %>%
  filter(value>0) %>%
  pivot_wider(names_from=name, values_from=value, values_fill=0) %>%
  left_join(all_meta %>% select(lib_id, genus, species, latitude, longitude)) %>%
  arrange(genus, species)

arbo_mat = arbo_table %>%
  select(any_of(arbovirus$virus_name)) %>%
  as.matrix()

rownames(arbo_mat) = arbo_table$species

RPM_col_fun = circlize::colorRamp2(c(0, ceiling(max(log10(arbo_mat+1)))), c("#F4F2ED", "#DD1332"))

Cairo::CairoPDF(file="output/Figure3b.pdf",       
                width=5,
                height=3.5)

plot(ComplexHeatmap::Heatmap(log10(t(arbo_mat)+1), 
                             column_order = 1:nrow(arbo_mat),
                             row_order = 1:ncol(arbo_mat),
                             col = RPM_col_fun,
                             name="log(RPM+1)"))
dev.off()


vvv = c("bunya0003")
arbo_table = select(rpm_table, lib_id, any_of(vvv)) %>% 
  pivot_longer(cols=2:ncol(.)) %>%
  filter(value>0) %>%
  pivot_wider(names_from=name, values_from=value, values_fill=0) %>%
  left_join(all_meta %>% select(lib_id, genus, species, latitude, longitude)) %>%
  arrange(genus, species)

map_data = arbo_table %>%
  pivot_longer(any_of(vvv)) %>%
  filter(value > 0)

china_map = sf::read_sf("map_data/Province.shp")

cluster_coord = sf::st_as_sf(
  map_data,
  coords = c("longitude", "latitude"),
  crs = 4326
)
f3a = ggplot() + 
  geom_sf(data=sf::st_transform(china_map, crs=4326)) + 
  geom_sf(aes(color=name), color="red", data=cluster_coord) +
  theme_bw() +
  ylab("Latitude") +
  xlab("Longitude") +
  ggsci::scale_fill_npg() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank())
f3a
ggsave("output/Figure3a.pdf", plot=f3a, width=6, height=5)
