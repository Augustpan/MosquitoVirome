library(tidyverse)

# load raw data
lib_meta <- read_csv("raw_data/lib_meta.csv")
lib_meta_o <- read_csv("output/lib_meta_merged.csv")
rpm_table <- read_csv("output/rpm_table_virus_masked.csv")
read_count <- read_tsv("raw_data/norRNA_reads.txt")
clim_data <- read_csv("raw_data/climate.csv")
lulc_data <- read_csv("raw_data/lulc_data.csv")
eco_data <- read_csv("raw_data/China.eco.csv")

vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))

# reordering datasets
rpm_table <- arrange(rpm_table, lib_id)
all_meta <- lib_meta %>%
  left_join(read_count) %>%
  left_join(clim_data) %>%
  left_join(lulc_data) %>%
  left_join(eco_data) %>%
  arrange(lib_id)

geo_cluster <- SoDA::geoXY(
  all_meta$latitude,
  all_meta$longitude,
  unit = 1000) %>%
  dist() %>%
  hclust(method = "average") %>%
  cutree(h=100)

# assign cluster info to each sample
all_meta$geo_cluster <- paste0("c", geo_cluster)

#### per individual analysis ####

# calculate virus richness **per individual**
richness <- rpm_table %>%
  select(-lib_id) %>%
  apply(1,function(x) {sum(x>0)})

all_meta$richness_per_ind <- richness

clim_nona = select(all_meta, lib_id, starts_with("bioclim_var")) %>%
  na.omit()
clim_pc = clim_nona %>%
  select(-lib_id) %>%
  scale() %>%
  prcomp()
clim_nona$CPC1 =  clim_pc$x[,1]
clim_nona$CPC2 =  clim_pc$x[,2]
clim_nona$CPC3 =  clim_pc$x[,3]

all_nona = all_meta %>% select(-starts_with("bioclim")) %>% right_join(clim_nona)


df = rpm_table %>%
  pivot_longer(cols=2:ncol(.), names_to="virus_name", values_to="rpm") %>%
  filter(rpm > 0) %>%
  left_join(select(lib_meta_o, lib_id, geo_cluster, species), by="lib_id") %>%
  group_by(geo_cluster, species, virus_name) %>%
  summarise(rpm = sum(rpm)) %>%
  ungroup() %>%
  pivot_wider(names_from="virus_name", values_from="rpm", values_fill=0)

comm_mat = select(df, geo_cluster)

climate_vars <- c("CPC1", "CPC2", "CPC3")
land_use_vars <- c(
  "Gridded_Mammal_Richness_2015",
  "log10(population_countKm2_2020)",
  "anthrome_PC1",
  "anthrome_PC2"
)
host_vars <- c("species")

candidate_vars <- c(climate_vars, land_use_vars, host_vars)
all_nona$`log10(population_countKm2_2020)` = log10(all_nona$population_countKm2_2020)
meta = all_nona %>%
  select(geo_cluster, all_of(candidate_vars), -species) %>%
  group_by(geo_cluster) %>%
  summarise_all(~mean(.x, na.rm=T)) %>%
  ungroup()

ll = all_nona %>%
  select(geo_cluster, all_of(candidate_vars)) %>%
  group_by(geo_cluster, species) %>%
  summarise() %>% ungroup() %>%
  left_join(meta)
meta = ll
rs = select(meta, geo_cluster, species) %>%
  left_join(df) %>%
  select(3:ncol(.)) %>%
  is.na() %>%
  rowSums()

ndf = select(meta, geo_cluster, species) %>%
  left_join(df)
ndf = ndf[rs<392,]
comm = select(ndf, -1, -2)
nmeta = meta[rs<392,]

nmeta = na.omit(nmeta)

ndf = select(nmeta, geo_cluster, species) %>%
  left_join(ndf)

comm = select(ndf, -1, -2)
mod0 = vegan::dbrda(comm ~ 1, data=nmeta, na.action = na.omit, dist="bray")
mod1 = vegan::dbrda(comm ~ ., data=nmeta, na.action = na.omit, dist="bray")
mods = vegan::ordistep(mod0, 
                         scope = formula(mod1), 
                         direction = "forward")

