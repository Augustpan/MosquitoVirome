library(tidyverse)
library(vegan)

# load raw data
lib_meta = read_csv("raw_data/lib_meta.csv")
rpm_table = read_csv("output/rpm_table_virus_masked.csv")
read_count = read_tsv("raw_data/norRNA_reads.txt")
clim_data = read_csv("raw_data/climate.csv")
#cst_profile = read_csv("output/cst_hiconf.csv")

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

tsne = Rtsne::Rtsne(rpm_table %>% select(-lib_id) %>% log1p(), 
                    check_duplicates = F,
                    max_iter = 5000,
                    perplexity=40)

all_meta$genus_alt = ifelse(all_meta$genus %in% c("Lutzia","Mansonia","Coquillettidia","Mimomyia"), 
                            "Others", all_meta$genus)
color_mapping = setNames(
  c("#3F5388", "#4DBBD6", "#E64B37", "#019F87", "#F29B7F", "#F29B7F", "#F29B7F", "#F29B7F", "#F29B7F"),
  c("Armigeres", "Anopheles", "Culex", "Aedes", "Mansonia", "Mimomyia", "Coquillettidia", "Lutzia", "Others")
)
f4a = ggplot(aes(x=tsne$Y[,1], y=tsne$Y[,2], color=all_meta$genus_alt), data=all_meta) +
  geom_point(size=1) + 
  scale_color_manual(values=color_mapping) +
  #ggsci::scale_color_npg() + 
  theme_bw() +
  xlab("T-SNE axis1") +
  ylab("T-SNE axis2")

ggsave("output/Figure4a.pdf", plot=f4a, width=5, height=3)
