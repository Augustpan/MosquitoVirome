library(tidyverse)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

rpm_table = read_csv("output/rpm_table_virus_masked.csv")
vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))

vtable = rpm_table %>%
  select_if(function(x) {
    ifelse(is.character(x), T, sum(x) > 0)
  }) # drop empty columns

virus_name_arr = colnames(select(vtable, -lib_id))
lib_id_arr = vtable$lib_id

vmeta = read_csv("output/virus_meta_merged.csv")
vmeta = filter(vmeta, is_core)
lib_meta = read_csv("raw_data/lib_meta.csv") %>%
  arrange(genus, species, city, year) %>%
  mutate(display_order = 1:nrow(.))

vmat = t(select(vtable, all_of(virus_name_arr)))
colnames(vmat) = vtable$lib_id

# reordering
vmeta_ord = arrange(vmeta, display_order)
meta_ord = arrange(lib_meta, display_order)
vmat_ord = log10(vmat[vmeta_ord$virus_name, meta_ord$lib_id] + 1)

# set color mapping
genus = sort(lib_meta$genus %>% unique())
genus_color = brewer.pal(length(genus), "Set1")

family = unique(vmeta_ord$virus_group)
family_color = c(brewer.pal(12, "Set3"),
                 brewer.pal(8, "Set2"),
                 brewer.pal(5, "Set1"))

RPM_col_fun = colorRamp2(c(0, ceiling(max(vmat_ord))), c("#F4F2ED", "#DD1332"))

# heatmap annotations
col_ano = columnAnnotation(
  host_genus = anno_simple(
    meta_ord$genus,
    col = setNames(genus_color, genus),
    height = unit(0.3, "cm")
  ),
  gap = unit(1, "mm")
)
row_ano = rowAnnotation(virus_family = anno_simple(
  vmeta_ord$virus_group,
  col = setNames(family_color, family),
  width = unit(0.3, "cm")
))

# draw heatmap LONG
pdf(file = "output/virus_heatmap_long.pdf",
    width = 12,
    height = 24)
ht = Heatmap(
  vmat_ord,
  cluster_columns = F,
  cluster_rows = F,
  row_order = 1:nrow(vmat_ord),
  column_order = 1:ncol(vmat_ord),
  column_split = factor(meta_ord$genus, levels = genus, ordered = T),
  name="log(RPM+1)",
  row_split = factor(vmeta_ord$virus_group, levels = family, ordered = T),
  col = RPM_col_fun,
  row_title_rot = 0,
  column_title_rot = 30,
  row_names_gp = gpar(fontsize = 6),
  top_annotation = col_ano,
  left_annotation = row_ano,
  column_labels = rep("", ncol(vmat_ord)),
  use_raster = F
)

l1 = Legend(
  labels = genus,
  title = "Host genus",
  legend_gp = gpar(fill = genus_color)
)
l4 = Legend(
  labels = family,
  title = "Virus taxa",
  legend_gp = gpar(fill = family_color)
)
draw(ht, annotation_legend_list = list(l1, l4))
dev.off()

# draw heatmap SHORT
pdf(file = "output/virus_heatmap.pdf",
    width = 12,
    height = 10)
ht = Heatmap(
  vmat_ord,
  cluster_columns = F,
  cluster_rows = F,
  row_order = 1:nrow(vmat_ord),
  column_order = 1:ncol(vmat_ord),
  column_split = factor(meta_ord$genus, levels = genus, ordered = T),
  row_split = factor(vmeta_ord$virus_group, levels = family, ordered = T),
  name="log(RPM+1)",
  col = RPM_col_fun,
  row_title_rot = 0,
  column_title_rot = 30,
  row_names_gp = gpar(fontsize = 6),
  top_annotation = col_ano,
  left_annotation = row_ano,
  row_labels = rep("", nrow(vmat_ord)),
  column_labels = rep("", ncol(vmat_ord)),
  use_raster = F
)

l1 = Legend(
  labels = genus,
  title = "Host genus",
  legend_gp = gpar(fill = genus_color)
)
l4 = Legend(
  labels = family,
  title = "Virus taxa",
  legend_gp = gpar(fill = family_color)
)
draw(ht, annotation_legend_list = list(l1, l4))
dev.off()
