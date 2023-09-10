library(tidyverse)
library(igraph)


library(tidyverse)

rpm_table = read_csv("output/rpm_table_virus_masked.csv") 
vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome)) %>%
  mutate(across(2:ncol(.), ~ifelse(.x>1, .x, 0)))

lib_meta = read_csv("raw_data/lib_meta.csv")

virus_names = colnames(rpm_table)[2:ncol(rpm_table)]
host_df = select(lib_meta, lib_id, genus, species)

incidence_matrix = rpm_table %>%
  left_join(host_df, by="lib_id") %>%
  relocate("species", .after=1) %>% 
  relocate("genus", .after=1) %>%
  mutate(across(all_of(virus_names), ~ifelse(.x>0, 1, 0)))

num_host_species = incidence_matrix %>%
  select(-lib_id, -genus) %>%
  group_by(species) %>%
  summarise_all(max) %>%
  select(-species) %>%
  colSums()

num_host_genus = incidence_matrix %>%
  select(-lib_id, -species) %>%
  group_by(genus) %>%
  summarise_all(max) %>%
  select(-genus) %>%
  colSums()

assertthat::are_equal(names(num_host_genus), names(num_host_species))
virus_host_stat = tibble(
  virus_name = names(num_host_species),
  num_host_species,
  num_host_genus
)

drop_empty = function(data) {
  cf = colSums(data) > 0
  rf = rowSums(data[,cf]) > 0
  return(list(data=data[rf,cf], row_filter=rf, col_filter=cf))
}

lib_meta = read_csv("raw_data/lib_meta.csv")
virus_name_arr = colnames(select(rpm_table, -lib_id))
virus_meta =  tibble(virus_name = virus_name_arr) %>%
  arrange(virus_name_arr) %>%
  mutate(order=1:nrow(.)) %>%
  mutate(virus_family = str_split(virus_name, "[\\d+_]", simplify = T)[,1])
vtable = select(read_csv("output/rpm_table_virus_masked.csv") , lib_id, all_of(core_virome)) %>%
  select_if(function(x) { ifelse(is.character(x), T, sum(x)>0)}) # drop empty columns
rpm_info = vtable

n_sample = lib_meta %>%
  count(species) %>% 
  arrange(desc(n))

top10_sp = n_sample$species[1:30]

stat_by_sp = lib_meta %>%
  select(lib_id, genus=genus, species=species) %>%
  right_join(rpm_info) %>%
  select(-lib_id) %>%
  group_by(genus, species) %>%
  summarise_all(~sum(.>0)) %>%
  ungroup()

mat = as.matrix(select(stat_by_sp, -genus, -species))
rownames(mat) = stat_by_sp$species
mat = mat[top10_sp,]
mat = mat[,colSums(mat)>0]
mat = mat[rowSums(mat)>0,]
vhg = graph.incidence(mat, weighted =T)
virus_meta = virus_host_stat %>% 
  mutate(cls = ifelse(num_host_species==1, "sps", ifelse(num_host_genus>1, "mpg", "ges"))) %>%
  right_join(virus_meta)
tmp = left_join(tibble(virus_name = colnames(select(stat_by_sp, -species, -genus))), virus_meta) %>%
  mutate(col = ifelse(cls=="sps", "#00FF00AA", ifelse(cls=="ges", "#FFFF00DD", "#FF0000DD")))

col_map = setNames(c(tmp$col, rep("#0000FF", length(stat_by_sp$species))), c(tmp$virus_name, stat_by_sp$species))
sz_map = setNames(c(rep(3, nrow(tmp)), rep(6, length(stat_by_sp$species))),  
                  c(tmp$virus_name, stat_by_sp$species))
dgv = mat %>% colSums()
sz = (rank(dgv, ties.method = "min") / max(rank(dgv, ties.method = "min")))#*0.5 + 1
alpha_map = setNames(
  c(sz[tmp$virus_name]*4, rep(3, length(stat_by_sp$species))),  
  c(tmp$virus_name, stat_by_sp$species)
)

flt = V(vhg)$name %in% stat_by_sp$species


lab = setNames(c(V(vhg)$name[flt], rep("", length(V(vhg)$name))[!flt]), c(V(vhg)$name[flt], V(vhg)$name[!flt]))




#col_fun = circlize::colorRamp2(c(1, max(degree(vhg)[tmp$virus_name], na.rm=T)), c("#aa00aa", "#FF0000"))
#col_map = setNames(c(col_fun(degree(vhg)[tmp$virus_name]), rep("#FFFF00", length(stat_by_sp$species))),
#                   c(tmp$virus_name, stat_by_sp$species))


pdf(file="output/Figure4b.pdf", width=10, height=10)
plot(vhg, 
     vertex.size=alpha_map[V(vhg)$name]*1.2, 
     vertex.color=col_map[V(vhg)$name], 
     vertex.label="",#lab,
     edge.color = "#00000020",
     edge.width=0.5)
dev.off()

