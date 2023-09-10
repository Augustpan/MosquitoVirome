library(tidyverse)

rpm_table = read_csv("output/rpm_table_virus_masked.csv")
lib_meta = read_csv("output/lib_meta_merged.csv")
virus_meta = read_csv("output/virus_meta_merged.csv")
mantel_corr = read_csv("output/mantel_corr.csv")

median_dist = function(latitude, longitude) {
  
  distance_mat = SoDA::geoXY(latitude, 
                             longitude, 
                             unit=1000) %>%
    dist()
  
  median_dist = median(distance_mat)
  median_dist
}

x = rpm_table %>%
  pivot_longer(cols = 2:ncol(.), 
               names_to = "virus_name", 
               values_to = "rpm") %>%
  filter(rpm > 0) %>%
  left_join(select(lib_meta, lib_id, latitude, longitude, genus, species)) %>%
  group_by(virus_name) %>%
  summarise(n=n(), 
            med_dist=median_dist(latitude, longitude), 
            n_species=length(unique(species)),
            n_genus=length(unique(genus))) %>%
  full_join(virus_meta) %>%
  full_join(mantel_corr) %>%
  mutate(is_endemic = ifelse(med_dist < 1000 | is.na(med_dist), "endemic", "wide-spread"),
         `co-divergence` = ifelse(p < 0.05, "co-diverge", NA),
         `distance decay` = ifelse(ps< 0.05, "yes", NA)) %>%
  select(
    virus_name,
    virus_group,
    baltimore_class,
    length,
    is_core,
    uncertainty,
    n, n_species,n_genus,`spatial coverage`=med_dist,
    `Mantel r (host)` = r,
    `P value (host)` = p,
    `Mantel r (space)` = rs,
    `P value (space)` = ps,
    is_endemic,`co-divergence`,`distance decay`
    
  )

host_rank = rpm_table %>%
  select(lib_id, all_of(filter(virus_meta, is_core)$virus_name)) %>%
  pivot_longer(cols = 2:ncol(.), 
               names_to = "virus_name", 
               values_to = "rpm") %>%
  filter(rpm > 0) %>%
  left_join(select(lib_meta, lib_id, species, genus)) %>%
  group_by(virus_name, species) %>%
  summarise(max_rpm = max(rpm), npos = n()) %>%
  group_by(virus_name) %>%
  group_modify(~tibble(host_rank_rpm=order(.x$max_rpm, decreasing = T), 
                       host_rank_pos=order(.x$npos, decreasing = T), 
                       species=.x$species, 
                       max_rpm=.x$max_rpm,
                       npos=.x$npos))


df = rpm_table %>%
  select(lib_id, all_of(filter(virus_meta, is_core)$virus_name)) %>%
  pivot_longer(cols = 2:ncol(.), 
               names_to = "virus_name", 
               values_to = "rpm") %>%
  filter(rpm > 0) %>%
  left_join(select(lib_meta, lib_id, species, genus)) %>%
  group_by(virus_name) %>%
  summarise(nrpm1k = length(unique(species[rpm > 1000])),
            median_rpm = median(rpm))

qplot(df$nrpm1k, df$max_rpm)
ggplot(host_rank, aes(x=log10(max_rpm))) +
  geom_histogram() +
  facet_wrap(~host_rank_pos)

ggplot(host_rank, aes(x=as.factor(host_rank_pos), y=log10(max_rpm))) +
  geom_point(size=0.1, position = flexplot::position_jitterd(width=0.4, height=0))

ggplot(host_rank, aes(x=as.factor(host_rank_rpm), y=log10(max_rpm))) +
  geom_point(size=0.1, position = flexplot::position_jitterd(width=0.4, height=0))

write_excel_csv(x, "x2.csv", na="")


s = rpm_table %>%
  pivot_longer(cols = 2:ncol(.), 
               names_to = "virus_name", 
               values_to = "rpm") %>%
  filter(rpm > 0) %>%
  left_join(select(virus_meta, virus_name, is_core), 
            by="virus_name") %>%
  group_by(lib_id) %>%
  summarise(ncore = sum(is_core),
            ntotal = sum(rpm>0)) 
  