library(tidyverse)
library(vegan)
library(SoDA)
library(ecodist)
library(ape)

drop_empty = function(data) {
  cf = colSums(data) > 0
  rf = rowSums(data[,cf]) > 0
  return(list(data=data[rf,cf], row_filter=rf, col_filter=cf))
}

lib_meta = read_csv("raw_data/lib_meta.csv")
clim = read_csv("raw_data/climate.csv")
lib_meta = left_join(lib_meta, clim)

rpm_table = read_csv("output/rpm_table_virus_masked.csv")
vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))

vtable = rpm_table %>%
  select_if(function(x) { ifelse(is.character(x), T, sum(x)>0)}) %>%
  mutate(across(2:ncol(.), ~ifelse(.x>1, .x, 0)))
virus_name_arr = colnames(select(vtable, -lib_id))
virus_meta =  tibble(virus_name = virus_name_arr) %>%
  arrange(virus_name_arr) %>%
  mutate(order=1:nrow(.)) %>%
  mutate(virus_family = str_split(virus_name, "[\\d+_]", simplify = T)[,1])
rpm_info = vtable

coi_tree = read.tree("raw_data/non_contaminated_COI.trimal.aln.treefile")

mat = select(rpm_info, -lib_id)

ret = drop_empty(mat)


mat_clean = ret$data

lib_id_names = rpm_info$lib_id[ret$row_filter]
virus_names = colnames(mat_clean)

lib_id_has_coi = intersect(coi_tree$tip.label, lib_id_names)
flt = lib_id_names %in% lib_id_has_coi
mat_has_coi = mat_clean[flt,]
mati = drop_empty(mat_has_coi)$data
metai = tibble(lib_id = lib_id_names) %>%
  left_join(lib_meta)
metai = metai[flt,]

virus_dist_bc = vegdist(mati, "bray")
virus_dist_jb = vegdist(mati, "jaccard", binary = T)

xy = geoXY(metai$latitude, metai$longitude, unit=1000)
rownames(xy) = metai$lib_id
geo_dist = dist(xy)

coi_dist = cophenetic(coi_tree)
coi_dist = as.dist(coi_dist[metai$lib_id,metai$lib_id])

cophenetic(coi_tree) %>%
  reshape2::melt(cophenetic(coi_tree)[metai$lib_id,metai$lib_id], 
                 varnames = c("lib1", "lib2"), 
                 value.name = "coi_dist") %>%
  as_tibble() %>%
  left_join(select(metai, lib1=lib_id, s1=species, g1=genus)) %>%
  left_join(select(metai, lib2=lib_id, s2=species, g2=genus)) %>%
  na.omit() %>%
  mutate(ss = s1 == s2, sg = g1 == g2) %>%
  group_by(ss) %>%
  summarise(mean_dist = mean(coi_dist),
            median_dist = median(coi_dist),
            q75_dist = quantile(coi_dist, 0.75),
            q25_dist = quantile(coi_dist, 0.25))

time_dist = dist(metai$year)

clim_dist = dist(select(metai, starts_with("bioclim_var")))

tmp_data = data.frame(
  virus_dist_jb = as.vector(virus_dist_jb),
  virus_dist_bc = as.vector(virus_dist_bc),
  coi_dist = as.vector(coi_dist),
  geo_dist = as.vector(geo_dist),
  time_dist = as.vector(time_dist),
  clim_dist = as.vector(clim_dist)
)

if (F) {
  a = mantel(virus_dist_jb ~ coi_dist + geo_dist + time_dist, data=tmp_data, mrank=T)
  b = mantel(virus_dist_jb ~ geo_dist + coi_dist + time_dist, data=tmp_data, mrank=T)
  c = mantel(virus_dist_jb ~ time_dist + coi_dist + geo_dist, data=tmp_data, mrank=T)
  
  d = mantel(virus_dist_bc ~ coi_dist + geo_dist + time_dist, data=tmp_data, mrank=T)
  e = mantel(virus_dist_bc ~ geo_dist + coi_dist + time_dist, data=tmp_data, mrank=T)
  f = mantel(virus_dist_bc ~ time_dist + coi_dist + geo_dist, data=tmp_data, mrank=T)
  
  mantel_results = rbind(
    jb_v_coi = a,
    jb_v_geo = b,
    jb_v_tim = c,
    bc_v_coi = d,
    bc_v_geo = e,
    bc_v_tim = f
  )
  write.csv(mantel_results, "mantel_results.csv")
  
}

tbl_host_sp = rpm_info %>%
  left_join(select(lib_meta, lib_id, species), by="lib_id") %>%
  pivot_longer(cols=all_of(colnames(rpm_info)[2:ncol(rpm_info)]), names_to="virus_name", values_to="RPM") %>%
  filter(RPM > 0) %>%
  distinct(lib_id, virus_name) %>%
  table() 
vs_host_sp_df = proxy::dist(tbl_host_sp, function(x, y) {sum((x>0)&(y>0))}) %>%
  as.matrix() %>%
  reshape2::melt(varnames = c("host_sp1", "host_sp2"), value.name = "shared_virus")

cdm = as.matrix(clim_dist)
rownames(cdm) = metai$lib_id
colnames(cdm) = metai$lib_id
clim_dist_df = reshape2::melt(cdm, varnames = c("host_sp1", "host_sp2"), value.name = "clim_dist")
coi_dist_df = reshape2::melt(as.matrix(coi_dist), varnames = c("host_sp1", "host_sp2"), value.name = "coi_dist")
geo_dist_df = reshape2::melt(as.matrix(geo_dist), varnames = c("host_sp1", "host_sp2"), value.name = "geo_dist")
#vbc_dist_df = reshape2::melt(as.matrix(virus_dist_bc), varnames = c("host_sp1", "host_sp2"), value.name = "vbc_dist")
#vjb_dist_df = reshape2::melt(as.matrix(virus_dist_jb), varnames = c("host_sp1", "host_sp2"), value.name = "vjb_dist")
all_dist_df = vs_host_sp_df %>%
  left_join(coi_dist_df) %>%
  left_join(geo_dist_df) %>%
  left_join(clim_dist_df) %>%
  #    left_join(vjb_dist_df) %>%
  #    left_join(vbc_dist_df) %>%
  na.omit() %>%
  filter(host_sp1 != host_sp2) %>%
  filter(coi_dist > 0)

M0 = glm(shared_virus ~ coi_dist + geo_dist, data=all_dist_df, family="poisson")

np = 5001
host_x = (0:(np-1))/(np-1)
geo_x = (0:(np-1))/(np-1) *5000
clim_x = host_x
host_pred = predict(M0, type="response", se.fit = T, newdata = tibble(coi_dist=host_x, geo_dist=rep(0,np)))
geo_pred =  predict(M0, type="response", se.fit=T, newdata = tibble(coi_dist=rep(0,np), geo_dist=geo_x))

sampled_data = sample_n(all_dist_df, 5000)
x_h = c(host_x,rev(host_x))
y_h = c(host_pred$fit + 2*host_pred$se.fit, rev(host_pred$fit-2*host_pred$se.fit))
g1 = ggplot() + 
  geom_jitter(aes(x=coi_dist, y=shared_virus), 
              data=sampled_data, 
              height=0.1, alpha=0.9, color="gray", shape=1) %>%
  ggrastr::rasterize(dpi=300) +
  geom_polygon(aes(x_h,y_h), fill="blue", alpha=0.5) +
  geom_line(aes(x=host_x, y=host_pred$fit)) +
  xlab("Phylogenetic distance") +
  ylab("Shared num. of virus") +
  theme_bw()

x_g = c(geo_x,rev(geo_x))
y_g = c(geo_pred$fit + 2*geo_pred$se.fit, rev(geo_pred$fit-2*geo_pred$se.fit))
g2 = ggplot() + 
  geom_jitter(aes(x=geo_dist, y=shared_virus), 
              data=sampled_data, 
              height=0.1, alpha=0.9, color="gray", shape=1) %>%
  ggrastr::rasterize(dpi=300) +
  geom_polygon(aes(x_g,y_g), fill="blue", alpha=0.5) +
  geom_line(aes(x=geo_x, y=geo_pred$fit)) +
  xlab("Geographic distance (km)") +
  ylab("Shared num. of virus") +
  theme_bw()

f4d = g1
f4e = g2
cowplot::plot_grid(g1, g2, nrow=2)
ggsave("output/Figure4de.pdf", width=3, height=4.6)
