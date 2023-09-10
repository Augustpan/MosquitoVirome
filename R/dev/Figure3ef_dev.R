library(tidyverse)

lib_meta = read_csv("output/lib_meta_merged.csv")
rpm_table = read_csv("output/rpm_table_virus_masked.csv")

vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))

geo_cluster = SoDA::geoXY(lib_meta$latitude, lib_meta$longitude, unit = 1000) %>%
  dist() %>%
  hclust(method = "average") %>%
  cutree(h=100)

lib_meta$geo_cluster = paste0("c", geo_cluster)

read_count = read_tsv("raw_data/norRNA_reads.txt")
richness = rpm_table %>%
  select(-lib_id) %>%
  apply(1,function(x) {sum(x>0)})
vload = rpm_table %>%
  select(-lib_id) %>%
  apply(1,function(x) {sum(x)})

x = tibble(lib_id = rpm_table$lib_id, richness=richness, vload=vload) %>%
  left_join(lib_meta) %>%
  left_join(read_count)

sp_rank = lib_meta %>% count(species) %>% arrange(desc(n))

fit = glm(richness ~ species + dbMEM.1 + dbMEM.2 + dbMEM.3 + dbMEM.4 +dbMEM.5 + dbMEM.6 + dbMEM.7 + dbMEM.8 +dbMEM.9 + dbMEM.10 + dbMEM.11 + dbMEM.12 +dbMEM.13 + dbMEM.14 + dbMEM.15 + dbMEM.16 + dbMEM.17 + year + log(norRNA_reads), family="poisson", data=filter(x, species %in% sp_rank$species[1:10]))
drop1(fit, test="Chisq")

fit = lm(log(vload+1) ~ species + dbMEM.1 + dbMEM.2 + dbMEM.3 + dbMEM.4 +dbMEM.5 + dbMEM.6 + dbMEM.7 + dbMEM.8 +dbMEM.9 + dbMEM.10 + dbMEM.11 + dbMEM.12 +dbMEM.13 + dbMEM.14 + dbMEM.15 + dbMEM.16 + dbMEM.17 + year + log(norRNA_reads), data=filter(x, species %in% sp_rank$species[1:10]))
anova(fit)

fit = glm(richness ~ species + dbMEM.2 + dbMEM.3 + dbMEM.4 +dbMEM.5 + dbMEM.6 + dbMEM.7 + dbMEM.8 + dbMEM.11 + dbMEM.12 +dbMEM.13 + dbMEM.14 + dbMEM.16 + year + log(norRNA_reads), family="poisson", data=filter(x, species %in% sp_rank$species[1:10]))
fit2 = lm(log(vload+1) ~ species + dbMEM.1 + dbMEM.2 + dbMEM.3 + dbMEM.4 +dbMEM.5 + dbMEM.6 + dbMEM.7 + dbMEM.8 + dbMEM.12 +dbMEM.13 + year + log(norRNA_reads), data=filter(x, species %in% sp_rank$species[1:10]))

species_fitted = predict(fit, terms = "species",type = "terms", se=T)
pred_df = tibble(lib_id =filter(x, species %in% sp_rank$species[1:10])$lib_id, 
                 species=filter(x, species %in% sp_rank$species[1:10])$species, 
                 sp_effect = species_fitted$fit[,1],
                 se = species_fitted$se.fit[,1]) %>%
  group_by(species) %>%
  summarise(sp_effect = mean(sp_effect),
            se = mean(se))

f3e = ggplot(aes(x=species, y=sp_effect), data=pred_df) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(number of virus per individual)")

ggsave("output/Figure3e.pdf", plot=f3e, width=5, height=4)


species_fitted = predict(fit2, terms = "species",type = "terms", se=T)
pred_df = tibble(lib_id =filter(x, species %in% sp_rank$species[1:10])$lib_id, 
                 species=filter(x, species %in% sp_rank$species[1:10])$species, 
                 sp_effect = species_fitted$fit[,1],
                 se = species_fitted$se.fit[,1]) %>%
  group_by(species) %>%
  summarise(sp_effect = mean(sp_effect),
            se = mean(se))

f3f = ggplot(aes(x=species, y=sp_effect), data=pred_df) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(viral RPM per individual)")
ggsave("output/Figure3f.pdf", plot=f3f, width=5, height=4)

# correlations
g1 = ggplot(aes(x=log1p(vload), y=richness), data=x) + 
  geom_point(size=1, shape=1) + 
  geom_smooth(method="lm") +
  xlab("log(total viral RPM+1)") +
  ylab("Num. of virus\nper individual") + 
  theme_bw() +
  ggpubr::stat_cor(method="spearman")
g2 =  ggplot(aes(x=log1p(norRNA_reads), y=richness), data=x) +
  geom_point(size=1, shape=1) +
  geom_smooth(method="lm") +
  xlab("log(total non-rRNA reads+1)") +
  ylab("Num. of virus\nper individual") + 
  theme_bw() +
  ggpubr::stat_cor(method="spearman")
g3 = ggplot(aes(x=log1p(norRNA_reads), y=log1p(vload)), data=x) + 
  geom_point(size=1, shape=1) +
  xlab("log(total non-rRNA reads+1)") +
  ylab("log(total viral RPM+1)") + 
  theme_bw() +
  ggpubr::stat_cor()

#### Figure S3 ####
cowplot::plot_grid(g1, g2, g3, nrow=2, align = "hv", labels=c("a", "b", "c"))
#ggsave("output/SI_depth_stat.pdf", width=4, height=7)
ggsave("output/SI_depth_stat.jpg",
       width=150, height=125, dpi=180, units="mm")