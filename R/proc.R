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

mantel_df = NULL
vdist_all = NULL
pmg_all = NULL
all_df = NULL
xx = c("partiti0042", "narna0001", "narna0005")
#for (v in top_v$virus_name) {
for (v in x) {
  
    tree_fname = str_glue("raw_data/viral_genomes_mapping/split_by_virus_p80/{v}.phy_phyml_tree.txt")
  

  if (file.exists(tree_fname) & file.size(tree_fname) > 0) {
    try({
        tree = tree_fname %>%
          ape::read.tree() %>%
          phytools::midpoint.root()

        vv =  cophenetic(tree)

          vn = str_extract(rownames(vv), "^.+_([0-9A-Z]+)_len\\d+_mapped\\d+$", group=1)
 
        

        if (length(vn) >= 1) {
          vdist = flatten_dist(vv, T) %>%
            mutate(virus_name = v,
                   Var1 = str_extract(Var1, "^.+_([0-9A-Z]+)_len\\d+_mapped\\d+$", group=1),
                   Var2 = str_extract(Var2, "^.+_([0-9A-Z]+)_len\\d+_mapped\\d+$", group=1))
          vdist_all = rbind(vdist_all, vdist)

          rownames(vv) = vn
          colnames(vv) = vn
          hh_v = hh[vn,vn]
          ss_v = ss[vn,vn]
          
          #pm1 = vegan::mantel.partial(as.dist(vv), as.dist(hh_v), as.dist(ss_v), method="spearman", parallel = 12)
          #pm2 = vegan::mantel.partial(as.dist(vv), as.dist(ss_v), as.dist(hh_v), method="spearman", parallel = 12)
          #pmg = ecodist::pmgram(as.dist(vv), as.dist(ss_v), as.dist(hh_v), stepsize=500)
          #pmg_all = rbind(pmg_all, mutate(pmg$mgram, virus=v, .before=1))
          #tdf = tibble(virus_name=v, species="all", n=length(vn), hcor = pm1$statistic, hp = pm1$signif, scor=pm2$statistic, sp = pm2$signif)
          #mantel_df = rbind(mantel_df, tdf)
          df = flatten_dist(vv, T) %>%
            rename(dist.v = dist) %>%
            left_join(flatten_dist(hh_v, T) %>% rename(dist.h = dist), by=c("Var1", "Var2")) %>%
            left_join(flatten_dist(ss_v, T) %>% rename(dist.s = dist), by=c("Var1", "Var2"))
          all_df = rbind(all_df, select(df, dist.v, dist.h, dist.s) %>% mutate(virus_name=v))
        } 
        # else {
          #tdf = tibble(virus_name=v, species="all", n=length(vn_sub), hcor = NA, hp = NA, scor = NA, sp = NA)
          #mantel_df = rbind(mantel_df, tdf)
        # }
        
        # for (sp in top_h$species[1:5] ) {
        #   sp_libs = filter(lib_meta, species == sp)$lib_id
        #   vn_sub = intersect(vn, sp_libs)
        #   if (length(vn_sub) >= 5) {
        #     vv_sub = vv[vn_sub, vn_sub]
        #     hh_sub = hh[vn_sub, vn_sub]
        #     ss_sub = ss[vn_sub, vn_sub]
        #     pms1 = vegan::mantel.partial(as.dist(vv_sub), as.dist(hh_sub), as.dist(ss_sub), method="spearman", parallel = 12)
        #     pms2 = vegan::mantel.partial(as.dist(vv_sub), as.dist(ss_sub), as.dist(hh_sub), method="spearman", parallel = 12)
        #     tdf = tibble(virus_name=v, species=sp, n=length(vn_sub), hcor = pms1$statistic, hp = pms1$signif, scor=pms2$statistic, sp = pms2$signif)
        #     mantel_df = rbind(mantel_df, tdf)
        #   } else {
        #     tdf = tibble(virus_name=v, species=sp, n=length(vn_sub), hcor = NA, hp = NA, scor = NA, sp = NA)
        #     mantel_df = rbind(mantel_df, tdf)
        #   }
        # }
      })
  }
}

# mantel_df %>%
#   group_by(virus_name) %>% 
#   summarise(host_codivergence=any(hp < 0.05, na.rm = T), 
#             spatial_divergence=any(sp < 0.05, na.rm=T)) %>%
#   left_join(vs, by="virus_name") %>%
#   write_csv("mantel_df_new.csv")

# all_dist_df = left_join(vdist_all, hs_dist, by=c("Var1", "Var2"))
# 
# subdf = all_dist_df %>%
#   filter(sp1 == sp2) %>%
#   filter(sp1 == top_h$species[1])

g1 = ggplot(data=filter(all_df, virus_name == "partiti0042"), aes(x=dist.s, y=dist.v)) + 
  ggrastr::rasterise(geom_point(size=0.5, alpha=0.1), dpi=300) + 
  geom_smooth(method="lm")
g2 = ggplot(data=filter(all_df, virus_name == "narna0001"), aes(x=dist.s, y=dist.v)) + 
  ggrastr::rasterise(geom_point(size=0.5, alpha=0.1), dpi=300) + 
  geom_smooth(method="lm")
g3 = ggplot(data=filter(all_df, virus_name == "narna0005"), aes(x=dist.s, y=dist.v)) + 
  ggrastr::rasterise(geom_point(size=0.5, alpha=0.1), dpi=300) + 
  geom_smooth(method="lm")

g = cowplot::plot_grid(g1,g2,g3, nrow=1)

ggsave(plot=g, filename="out.pdf", width=160, height=30, units="mm")
