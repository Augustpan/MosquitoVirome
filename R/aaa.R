flatten_dist = function(d) {
  dmat = as.matrix(d)
  ut = upper.tri(dmat)
  rownames(ut) = rownames(dmat)
  colnames(ut) = colnames(dmat)
  reshape2::melt(ut, value.name = "flag") %>%
    left_join(reshape2::melt(dmat, value.name = "dist")) %>%
    filter(flag) %>%
    select(-flag) %>%
    as_tibble()
}

lib_meta = read_csv("output/lib_meta_merged.csv")
cst_profile = read_csv("output/cst_vhconf.csv")

top_sp = arrange(count(lib_meta, species), desc(n))$species
mos = top_sp[5]
for (mos in top_sp[1:5]) {
  
  lib_meta_mos = filter(lib_meta, species %in% mos)
  
  geo_loc = select(lib_meta, geo_cluster, latitude, longitude) %>%
    group_by(geo_cluster) %>%
    summarise_all(mean)
  
  xy = SoDA::geoXY(geo_loc$latitude, geo_loc$longitude, unit=1000)
  rownames(xy) = geo_loc$geo_cluster
  sdist = flatten_dist(dist(xy))
  
  htree = ape::read.tree("raw_data/non_contaminated_COI.trimal.aln.treefile")
  hdist = flatten_dist(cophenetic(htree)) %>%
    left_join(select(lib_meta_mos, Var1=lib_id, c1=geo_cluster), by="Var1") %>% 
    left_join(select(lib_meta_mos, Var2=lib_id, c2=geo_cluster), by="Var2") %>% 
    na.omit() %>%
    
    group_by(c1, c2) %>% 
    summarise(min_hdist=min(dist), q50_hdist=median(dist), q10_hdist=quantile(dist,0.1)) %>%
    ungroup() %>%
    select(Var1=c1, Var2=c2, min_hdist, q50_hdist, q10_hdist)
  
  cstv = rpm_table %>%
    pivot_longer(cols=2:ncol(.), values_to="rpm", names_to="virus_name") %>%
    filter(rpm > 0) %>%
    filter(virus_name %in% cst_profile$virus_name)  %>%
    left_join(select(lib_meta_mos, lib_id, geo_cluster)) %>%
    na.omit() %>%
    group_by(virus_name, geo_cluster) %>%
    summarise(rpm = sum(rpm)) %>%
    pivot_wider(names_from="virus_name", values_from="rpm", values_fill=0)
  
  shared_virus_cst = select(cstv, -1) %>%
    vegan::designdist(method="a", abcd=T) %>%
    as.matrix()
  colnames(shared_virus_cst) = cstv$geo_cluster
  rownames(shared_virus_cst) = cstv$geo_cluster
  shared_virus_cst = flatten_dist(shared_virus_cst) %>%
    rename(shared_cst=dist)
  
  non_cstv =  rpm_table %>%
    pivot_longer(cols=2:ncol(.), values_to="rpm", names_to="virus_name") %>%
    filter(rpm > 0) %>%
    left_join(select(lib_meta_mos, lib_id, geo_cluster)) %>%
    na.omit() %>%
    group_by(virus_name, geo_cluster) %>%
    summarise(rpm = sum(rpm)) %>%
    pivot_wider(names_from="virus_name", values_from="rpm", values_fill=0)
  
  shared_virus_non_cst = select(non_cstv, -1) %>%
    vegan::designdist(method="a", abcd=T) %>%
    as.matrix()
  colnames(shared_virus_non_cst) = non_cstv$geo_cluster
  rownames(shared_virus_non_cst) = non_cstv$geo_cluster
  shared_virus_non_cst = flatten_dist(shared_virus_non_cst) %>%
    rename(shared_non_cst=dist)
  
  data_wide = sdist %>%
    left_join(hdist) %>%
    left_join(shared_virus_cst) %>%
    left_join(shared_virus_non_cst) %>%
    replace_na(list(shared_cst=0, shared_non_cst=0)) %>%
    na.omit()
  data = data_wide %>%
    pivot_longer(cols=shared_cst:shared_non_cst, names_to="cst", values_to="shared_virus")
  
  gh = ggplot(filter(data, cst=="shared_non_cst")) +
    geom_point(aes(x=min_hdist, y=shared_virus)) +
    geom_smooth(aes(x=min_hdist, y=shared_virus),
                method="glm", method.args=list(family="poisson"), se=F) +
    theme_bw() + 
    theme(legend.position = "none")
  
  gs = ggplot(filter(data, cst=="shared_non_cst")) +
    geom_point(aes(y=min_hdist, x=dist)) +
    geom_smooth(aes(y=min_hdist, x=dist), method="lm") + theme_bw()
  
  gvs = ggplot(filter(data, cst=="shared_non_cst")) +
    geom_point(aes(y=shared_virus, x=dist)) +
    geom_smooth(aes(y=shared_virus, x=dist), method="glm",  method.args=list(family="poisson"))+
    theme_bw() + 
    theme(legend.position = "none")
  
  g = cowplot::plot_grid(gh, gvs, gs, nrow=3)
  ggsave(plot=g, filename=str_glue("{mos}.pdf"), width=80, height=180, units="mm")
  
  
  data_with_loc = data_wide %>%
    left_join(rename(geo_loc, Var1=geo_cluster, lat1=latitude, lon1=longitude)) %>%
    left_join(rename(geo_loc, Var2=geo_cluster, lat2=latitude, lon2=longitude)) 
  power = function(x, y) x^y
  g1 = ggplot() + 
    geom_point(aes(x=longitude, y=latitude), data=geo_loc) +
    geom_segment(aes(x=lon1, xend=lon2, y=lat1, yend=lat2, 
                     alpha=hd, size=hd), 
                 data=mutate(data_with_loc, hd=ifelse(min_hdist<0.02, 1-min_hdist, 1-0.02))) +
    scale_alpha_continuous(range=c(0.1, 0.6)) +
    scale_size_continuous(range=c(0.01, 0.5)) + theme_bw()
  
  g2 = ggplot() + 
    geom_point(aes(x=longitude, y=latitude), data=geo_loc) +
    geom_segment(aes(x=lon1, xend=lon2, y=lat1, yend=lat2, 
                     alpha=shared_non_cst, size=shared_non_cst), 
                 data=mutate(data_with_loc)) +
    scale_alpha_continuous(range=c(0.1, 0.6)) +
    scale_size_continuous(range=c(0.01, 0.5)) + theme_bw()
  
  g3 = ggplot() + 
    geom_point(aes(x=longitude, y=latitude), data=geo_loc) +
    geom_segment(aes(x=lon1, xend=lon2, y=lat1, yend=lat2, 
                     alpha=shared_cst, size=shared_cst), 
                 data=mutate(data_with_loc)) +
    scale_alpha_continuous(range=c(0.1, 0.6)) +
    scale_size_continuous(range=c(0.01, 0.5)) + theme_bw()
  
  gg = cowplot::plot_grid(g1, g2, nrow=2, align="hv")
  ggsave(plot=gg, filename=str_glue("map_{mos}.pdf"), width=5.39, height=6.59)
}
