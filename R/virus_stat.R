library(tidyverse)

flatten_dist = function(d, tri=F) {
  dmat = as.matrix(d)
  if (tri) {
    ut = upper.tri(dmat)
  } else {
    ut = matrix(data=TRUE, nrow=nrow(dmat), ncol=ncol(dmat))
  }
  rownames(ut) = rownames(dmat)
  colnames(ut) = colnames(dmat)
  reshape2::melt(ut, value.name = "flag") %>%
    left_join(reshape2::melt(dmat, value.name = "dist"), by=c("Var1", "Var2")) %>%
    filter(flag) %>%
    select(-flag) %>%
    as_tibble()
}

lib_meta = read_csv("output/lib_meta_merged.csv")
rpm_table = read_csv("output/rpm_table_virus_masked.csv")

top_v = rpm_table %>%
  pivot_longer(cols=2:ncol(.), values_to="rpm", names_to="virus_name") %>%
  filter(rpm > 0) %>%
  count(virus_name) %>%
  arrange(desc(n))

top_h = count(lib_meta, species) %>% arrange(desc(n))
geo_loc = select(lib_meta, geo_cluster, latitude, longitude) %>%
  group_by(geo_cluster) %>%
  summarise_all(mean)

xy = SoDA::geoXY(lib_meta$latitude, lib_meta$longitude, unit=1000)
rownames(xy) = lib_meta$lib_id
ss = as.matrix(dist(xy))
sdist = flatten_dist(ss)

htree = ape::read.tree("raw_data/non_contaminated_COI.trimal.aln.treefile")
hh = cophenetic(htree)
hh = as.matrix(hh[lib_meta$lib_id, lib_meta$lib_id])
hdist = flatten_dist(hh)
hs_dist = full_join(sdist, hdist, by=c("Var1", "Var2")) %>%
  left_join(select(lib_meta, Var1=lib_id, sp1=species), by="Var1") %>%
  left_join(select(lib_meta, Var2=lib_id, sp2=species), by="Var2") %>%
  rename(sdist=dist.x, hdist=dist.y)

vs = rpm_table %>%
  pivot_longer(cols=2:ncol(.), values_to="rpm", names_to="virus_name") %>%
  filter(rpm > 0) %>%
  left_join(lib_meta, by="lib_id") %>%
  group_by(virus_name) %>%
  summarise(n_loc=length(unique(geo_cluster)), 
            n_ind=length(unique(lib_id)),
            n_cpip=length(unique(lib_id[species=="Culex pipiens"])),
            n_asub=length(unique(lib_id[species=="Armigeres subalbatus"])),
            n_aalb=length(unique(lib_id[species=="Aedes albopictus"])),
            n_ctri=length(unique(lib_id[species=="Culex tritaeniorhynchus"])),
            n_asin=length(unique(lib_id[species=="Anopheles sinensis"])),
            n_sp=length(unique(species)),
            n_ge=length(unique(genus)),
            min_lat=min(latitude), max_lat=max(latitude), lat_span_IQR=quantile(latitude, 0.75)-quantile(latitude, 0.25), 
            min_lon=min(longitude), max_long=max(longitude), lon_span_IQR=quantile(longitude, 0.75)-quantile(longitude, 0.25))