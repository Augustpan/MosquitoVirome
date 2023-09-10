library(tidyverse)

lib_meta = read_csv("output/lib_meta_merged.csv")
virus_meta = read_csv("output/virus_meta_merged.csv")
rpm_table = read_csv("output/rpm_table_virus_masked.csv")

core_virome = filter(virus_meta, is_core)$virus_name
rpm_core = select(rpm_table, lib_id, all_of(core_virome))

sp_count = lib_meta %>%
  count(species) %>%
  arrange(desc(n))

rarefy_mos = function(tbl, sample_size, rep=10) {
  f = function(r) {
    cs = select(sample_n(tbl, sample_size), -lib_id) %>%
      colSums()
    sum(cs>0)
  }
  s = map_dbl(1:rep, f)
  c(sample_size, mean(s), sd(s))
}
rarefy_mos_vec = Vectorize(rarefy_mos, "sample_size")
rarefy_mos_vec_wrap = function(tbl, sample_size, rep=10) {
  out = t(rarefy_mos_vec(tbl, sample_size)) %>%
    as.data.frame() %>%
    as_tibble()
  
  colnames(out) = c("sample_size", "mean", "sd")
  out
}
richness = rpm_core %>%
  left_join(select(lib_meta, lib_id, species), by="lib_id") %>%
  filter(species %in% sp_count$species[1:5]) %>%
  group_by(species) %>%
  group_modify(~rarefy_mos_vec_wrap(.x, seq(1, nrow(.x), 1)))

stat = richness %>% 
  group_by(species) %>% 
  summarise(x_max = max(sample_size), 
            y_max=max(mean)) %>%
  as.data.frame()

rownames(stat) = stat$species

g = ggplot(aes(x=sample_size, y=mean, color=species), data=richness) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                width=.1, color="blue", alpha=0.15) +
  geom_line(size=.3 )
  
for (sp in stat$species) {
  g = g + annotate(geom = "text", 
                   x =stat[sp,]$x_max, 
                   y = stat[sp,]$y_max,
                   label = sp,
                   size=2,
                   hjust=0)
}

f3g = g +
  xlab("Num. of virus species") +
  ylab("Num. of\nmosquito individuals") +
  theme_bw() + 
  theme(legend.position = "none")
