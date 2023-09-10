library(tidyverse)
library(vegan)
# load raw data
lib_meta <- read_csv("raw_data/lib_meta.csv")
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

#### per individual analysis ####d

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

all_nona = all_meta %>% select(-starts_with("bioclim")) %>% right_join(clim_nona) %>%
  mutate(log10_population_countKm2_2020 = log10(population_countKm2_2020))

climate_vars <- c("CPC1", "CPC2", "CPC3")
land_use_vars <- c(
  "Gridded_Mammal_Richness_2015",
  "log10_population_countKm2_2020",
  "anthrome_PC1",
  "anthrome_PC2"
)
host_vars <- c("species")

candidate_vars <- c(climate_vars, land_use_vars, host_vars)

top10 = (count(lib_meta, species) %>% arrange(desc(n)))$species[1:10]
envir = select(all_nona, all_of(candidate_vars), lib_id) %>% filter(species %in% top10) %>% na.omit()
comm = filter(rpm_table, lib_id %in% envir$lib_id) 

rs = select(comm, -lib_id) %>% rowSums()
to_drop = comm$lib_id[rs == 0]

comm = filter(comm, !(lib_id %in% to_drop)) %>% select(-lib_id)
envir = filter(envir, !(lib_id %in% to_drop))
envir$lib_id = NULL

comm_dist = designdist(comm, method = "a", abcd = TRUE)
dfunc = function(x, y) designdist(x, method = "a", abcd = TRUE)
mod0 = rda(log10(comm+1) ~ 1, data=envir)
mod1 = rda(log10(comm+1) ~ ., data=envir)
mods = ordiR2step(mod0, 
                  scope = formula(mod1), 
                  direction = "forward")

vp = varpart(log10(comm+1), 
        select(envir, CPC1, CPC2, CPC3), 
        select(envir, anthrome_PC1, anthrome_PC2, Gridded_Mammal_Richness_2015),
        select(envir, species))

envir = envir %>%
  mutate_if(is.double, scale)

mod0 = dbrda(comm ~ 1, data=envir, dist="jaccard")
mod1 = dbrda(comm ~ ., data=envir, dist="jaccard")
mods = ordiR2step(mod0, 
                  scope = formula(mod1), 
                  direction = "forward")

mods = dbrda(comm ~ species + CPC2 + Gridded_Mammal_Richness_2015 + anthrome_PC1 +
               CPC1 + CPC3 + anthrome_PC2, dist = "jaccard")

vp = varpart(vegdist(comm, "jaccard"), 
             select(envir, CPC1, CPC2, CPC3), 
             select(envir, anthrome_PC1, anthrome_PC2, Gridded_Mammal_Richness_2015),
             select(envir, species))

