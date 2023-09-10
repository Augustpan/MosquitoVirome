library(tidyverse)
library(ComplexHeatmap)
library(circlize)

rpm_table = read_csv("output/rpm_table_virus_masked.csv")
vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))

lib_meta = read_csv("raw_data/lib_meta.csv")
cst_profile = read_csv("output/cst_hiconf.csv")

tbl_rename = read_tsv("raw_data/renaming.txt")

cross_genera = cst_profile %>%
  distinct(genus, virus_name) %>% 
  count(virus_name) %>%
  filter(n>1)

cgt = cross_genera$virus_name

ddf = tibble(
  lid = rpm_table$lib_id,
  generalist = ((rpm_table %>%
                   select(all_of(cgt)))>0) %>%
    rowSums(),
  
  specialist = ((rpm_table %>%
                   select(-lib_id,-all_of(cgt)))>0) %>%
    rowSums()
)

virus_names = colnames(rpm_table)[2:ncol(rpm_table)]
host_df = select(lib_meta, lib_id, genus, species)

rpm_info_cst = rpm_table %>%
  left_join(host_df, by="lib_id") %>%
  relocate("species", .after=1) %>% 
  relocate("genus", .after=1) %>%
  pivot_longer(cols=4:ncol(.), names_to="virus_name", values_to="RPM") %>%
  filter(RPM > 0) %>%
  right_join(cst_profile) %>%
  pivot_wider(values_from ="RPM", values_fill=0, names_from="virus_name") 



sp_mat = rpm_info_cst %>%
  select(-genus, -lib_id) %>%
  group_by(species) %>%
  summarise_all(~ifelse(is.na(median(.x[.x!=0])), 0,median(.x[.x!=0])))


smat = as.matrix(select(sp_mat, -species) + 1) %>%
  log10()

rownames(smat)= sp_mat$species

smat = t(smat)
#RPM_col_fun = colorRamp2(c(0, ceiling(max(smat))), c("#F4F2ED", "#DD1332"))
RPM_col_fun = colorRamp2(c(0, 1), c("#F4F2ED", "#DD1332"))

pdf(file = "output/cst_heatmap_species.pdf",
    width = 12,
    height = 10)
plot(Heatmap(smat,         name = "median log(RPM+1)",col=RPM_col_fun))

dev.off()

ge_mat = rpm_info_cst %>%
  select(lib_id, genus, species, all_of(cgt)) %>%
  select(-species, -lib_id) %>%
  group_by(genus) %>%
  summarise_all(~ifelse(is.na(median(.x[.x!=0])), 0,median(.x[.x!=0])))

gmat = as.matrix(select(ge_mat, -genus) + 1) %>%
  log10()
rownames(gmat) = ge_mat$genus
gmat = t(gmat)

gca = vegan::cca(gmat)
gmat = gmat[order(gca$CA$u[,1]), order(gca$CA$v[,1])]

map_rename = setNames(tbl_rename$new_name, tbl_rename$virus_name)
rownames(gmat) = map_rename[rownames(gmat)]

RPM_col_fun = colorRamp2(c(0, ceiling(max(gmat))), c("#F4F2ED", "#DD1332"))

pdf(file = "output/cst_heatmap_genus.pdf",
    width = 5.5,
    height = 10)
Heatmap(
  gmat, 
  name = "median log(RPM+1)",
  col=RPM_col_fun,
  row_order = 1:nrow(gmat),
  column_order = 1:ncol(gmat),
) %>% plot()
dev.off()

# prevalence
ge_mat = rpm_info_cst %>%
  select(lib_id, genus, species, all_of(cgt)) %>%
  select(-species, -lib_id) %>%
  group_by(genus) %>%
  summarise_all(~ifelse(n()==0, 0, sum(.x>0)/n()))

gmat = as.matrix(select(ge_mat, -genus))
rownames(gmat)= ge_mat$genus

sp_mat= rpm_info_cst %>%
  select(-genus, -lib_id) %>%
  group_by(species) %>%
  summarise_all(~ifelse(n()==0, 0, sum(.x>0)/n()))
smat = as.matrix(select(sp_mat, -species))

rownames(smat)= sp_mat$species

smat = t(smat)
#RPM_col_fun = colorRamp2(c(0, ceiling(max(smat))), c("#F4F2ED", "#DD1332"))
RPM_col_fun = colorRamp2(c(0, 1), c("#F4F2ED", "#DD1332"))

pdf(file = "output/cst_heatmap_species_prevalence.pdf",
    width = 12,
    height = 10)
Heatmap(smat, 
        name = "Prevalence",
        col=RPM_col_fun) %>% plot()
dev.off()

gmat = t(gmat)
gmat[is.na(gmat)] = 0
gca = vegan::cca(gmat)
gmat = gmat[order(gca$CA$u[,1]), order(gca$CA$v[,1])]

RPM_col_fun = colorRamp2(c(0, ceiling(max(gmat))), c("#F4F2ED", "#DD1332"))

pdf(file = "output/cst_heatmap_genus_prevalence.pdf",
    width = 6,
    height = 10)
Heatmap(
  gmat, 
  name = "Prevalence",
  col=RPM_col_fun,
  row_order = 1:nrow(gmat),
  column_order = 1:ncol(gmat),
) %>% plot()
dev.off()
