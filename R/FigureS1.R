library(tidyverse)

lib_meta = read_csv("output/lib_meta_merged.csv")

geo_cluster = lib_meta$geo_cluster

cluster_centroids = lib_meta %>%
  select(geo_cluster, latitude, longitude) %>%
  group_by(geo_cluster) %>%
  summarise(lat=mean(latitude), lon=mean(longitude), nsample=n())

max_dist = function(tmp, ...) {
  md = SoDA::geoXY(tmp$latitude, tmp$longitude, unit=1000) %>%
    dist() %>% max()
  tibble(max_dist=md)
}

max_dist_within_clusters = lib_meta %>%
  select(geo_cluster, latitude, longitude) %>%
  group_by(geo_cluster) %>%
  group_modify(max_dist) %>%
  ungroup()

cluster_centroids = cluster_centroids %>%
  left_join(max_dist_within_clusters)

china_map = sf::read_sf("map_data/Province.shp")

cluster_coord = sf::st_as_sf(
  cluster_centroids,
  coords = c("lon", "lat"),
  crs = 4326
)

sp_count = lib_meta %>%
  count(species) %>% 
  arrange(desc(n)) %>%
  mutate(species = factor(species, levels=species))

draw_mos_distribution = function(data, group_var) {
  plot_coord = sf::st_as_sf(
    data,
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  ggplot() + 
    geom_sf(data=sf::st_transform(china_map, crs=4326)) + 
    geom_sf(data=plot_coord) +
    theme_bw() +
    ylab("Latitude") +
    xlab("Longitude") +
    ggsci::scale_fill_npg() +
    ggtitle(group_var[[1]]) +
    theme(plot.title = element_text(face="italic"))
}

#### Figure S2 ####
plot_data = lib_meta %>%
  filter(species %in% sp_count$species[1:5]) %>%
  select(species, geo_cluster) %>%
  left_join(cluster_centroids, by="geo_cluster") %>%
  group_by(species) %>%
  group_map(draw_mos_distribution)

combined = cowplot::plot_grid(plotlist = plot_data, nrow=2, ncol=3)
ggsave(filename="output/FigureS2.jpg", plot=combined, 
       width=220, height=120, units="mm", dpi=300)

#### Figure S1 ####
sp_count_plot = ggplot() +
  geom_bar(aes(x=species, y=n), data=sp_count[1:20,], stat="identity") +
  xlab("Top 20 most abundant mosquitos") + 
  ylab("Num. of samples") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, face="italic"))
ggsave(filename="output/FigureS1.jpg", plot=sp_count_plot, 
       width=100, height=80, units="mm", dpi=300)