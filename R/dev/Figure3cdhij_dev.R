library(tidyverse)

# load raw data
lib_meta <- read_csv("raw_data/lib_meta.csv")
rpm_table <- read_csv("output/rpm_table_virus_masked.csv")
read_count <- read_tsv("raw_data/norRNA_reads.txt")
clim_data <- read_csv("raw_data/climate.csv")
lulc_data <- read_csv("raw_data/lulc_data.csv")
eco_data <- read_csv("raw_data/China.eco.csv")

vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))

# reordering datasets
rpm_table <- arrange(rpm_table, lib_id)
all_meta <- lib_meta %>%
  left_join(read_count) %>%
  left_join(clim_data) %>%
  left_join(lulc_data) %>%
  left_join(eco_data) %>%
  arrange(lib_id)

geo_cluster <- SoDA::geoXY(
  all_meta$latitude,
  all_meta$longitude,
  unit = 1000) %>%
  dist() %>%
  hclust(method = "average") %>%
  cutree(h=100)

# assign cluster info to each sample
all_meta$geo_cluster <- paste0("c", geo_cluster)

#### per individual analysis ####

# calculate virus richness **per individual**
richness <- rpm_table %>%
  select(-lib_id) %>%
  apply(1,function(x) {sum(x>0)})

all_meta$richness_per_ind <- richness

clim_nona = select(all_meta, lib_id, starts_with("bioclim_var")) %>%
  na.omit()
clim_pc = clim_nona %>%
  select(-lib_id) %>%
  scale() %>%
  prcomp()
clim_nona$CPC1 =  clim_pc$x[,1]
clim_nona$CPC2 =  clim_pc$x[,2]
clim_nona$CPC3 =  clim_pc$x[,3]

all_nona = all_meta %>% select(-starts_with("bioclim")) %>% right_join(clim_nona)

f3h = ggplot(aes(y=richness_per_ind, x=bioclim_var1), data=all_nona) +
  geom_point(size=1, alpha=0.3) + 
  geom_smooth(method = "lm") +
  theme_bw()

f3i = ggplot(aes(y=richness_per_ind, x=bioclim_var12), data=all_nona) +
  geom_point(size=1, alpha=0.3) + 
  geom_smooth(method = "lm") +
  theme_bw()

f3j = ggplot(aes(y=richness_per_ind, x=Gridded_Mammal_Richness_2015), data=all_nona) +
  geom_point(size=1, alpha=0.3) + 
  geom_smooth(method = "lm") +
  theme_bw()

cowplot::plot_grid(f3h, f3i, f3j, nrow=1)
ggsave("output/Figure3hij.pdf", width=154, height=45, units="mm")

climate_vars <- c("CPC1", "CPC2", "CPC3")
land_use_vars <- c(
  "Gridded_Mammal_Richness_2015",
  "log10(population_countKm2_2020)",
  "anthrome_PC1",
  "anthrome_PC2"
)
host_vars <- c("species")

candidate_vars <- c(climate_vars, land_use_vars, host_vars)

all_nona <- all_nona %>%
  select(
    richness_per_ind,
    any_of(candidate_vars),
    population_countKm2_2020,
    year) %>%
  na.omit()

enum_formula <- function(n) {
  comb <- combinat::combn(candidate_vars, n)
  if (any(is.null(dim(comb)))) {
    ret <- paste0("richness_per_ind~", paste0(comb, collapse = "+"))
  } else {
    ret <- apply(comb, 2, function(x) {
      paste0(
        "richness_per_ind~as.factor(year)+",
        paste0(x, collapse = "+")
      )
    })
  }
  return(ret)
}

all_fitted <- purrr::map(1:length(candidate_vars),enum_formula) %>%
  reduce(c) %>%
  purrr::map(~glm(as.formula(.x), family="poisson", data=all_nona))

model_sel <- MuMIn::model.sel(all_fitted, rank=AIC)
model_index <- rownames(model_sel) %>% as.integer()

topN <- 10
var_explained <- all_fitted[model_index[1:topN]] %>%
  lapply(function(x) {
    drop_one <- drop1(x)
    drop_tbl <- drop_one %>%
      as_tibble() %>%
      mutate(term = rownames(drop_one), .before = 1) %>%
      mutate(dev_exp = (Deviance-x$deviance)/x$null.deviance) %>%
      select(-Df, -Deviance, -AIC)
    drop_tbl[1, 2] <- (x$null.deviance - x$deviance) / x$null.deviance
    drop_tbl
  }) %>%
  reduce(full_join, by="term")

colnames(var_explained) = c("term", paste0("dev_exp_", 1:topN))


f = paste0("richness_per_ind~", paste0(climate_vars, collapse="+")) %>%
  as.formula() %>%
  glm(family = "poisson", data=all_nona)
dev_explained_clim = (f$null.deviance-f$deviance) / f$null.deviance

f = paste0("richness_per_ind~", paste0(land_use_vars, collapse="+")) %>%
  as.formula() %>%
  glm(family = "poisson", data=all_nona)
dev_explained_land = (f$null.deviance-f$deviance) / f$null.deviance

f= paste0("richness_per_ind~", paste0(host_vars, collapse="+")) %>%
  as.formula() %>%
  glm(family = "poisson", data=all_nona)
dev_explained_host = (f$null.deviance-f$deviance) / f$null.deviance

var_mat = var_explained[2:10,2:11] %>% t()
colnames(var_mat) <- var_explained$term[2:10]

dev_prop = tibble(
  var = c("climate", "land use", "host species"),
  dev_exp = c(dev_explained_clim, dev_explained_land, dev_explained_host)
)
dev_prop = rbind(dev_prop, c("var"="unexplained", 
                             "dev_exp"=1-sum(dev_prop$dev_exp)))
dev_prop$dev_exp = as.double(dev_prop$dev_exp)


write.csv(dev_prop, "output/Figure3d_data.csv")

#### F3c ####

tbl_model_sel = as_tibble(var_mat) %>%
#  mutate(model_rank = 1:nrow(.)) %>%
  select(`Mosquito species` = species,
         `Climate PC1` = CPC1,
         `Climate PC2` = CPC2,
         `Climate PC3` = CPC3,
         `Land-use PC1` =  `anthrome_PC1`,
         `Land-use PC2` =  `anthrome_PC2`,
         `Mammal Richness` = Gridded_Mammal_Richness_2015,
         `log(population density)` = `log10(population_countKm2_2020)`,
         `Year of sample`=`as.factor(year)`) %>%
  mutate(deltaAIC=model_sel[1:topN,]$delta)

write_csv(tbl_model_sel, "output/Figure3c_data.csv")

#### F3d ####
plot_data = dev_prop %>%
  arrange(dev_exp) %>%
  mutate(prop = dev_exp / sum(dev_exp) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop)

plot_data$var = factor(plot_data$var, 
                       levels=c("climate", "land use", "host species", "unexplained"))

f3d = ggplot(plot_data) +
  geom_bar(aes(x="", y=-dev_exp, fill=var), color="black", stat="identity", width=1) +
  coord_polar("y", start= pi*0/180) +
  scale_fill_manual(values = c("climate"="#426db1", "land use"="#e47a3d", 
                               "host species"="#f4ba34", "unexplained"="#cac9c5")) + 
  theme(legend.position = "none", 
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
