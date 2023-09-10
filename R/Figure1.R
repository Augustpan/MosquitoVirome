library(tidyverse)

draw_figure1a = function(lib_meta, china_map) {
  max_dist = function(tmp, ...) {
    md = SoDA::geoXY(tmp$latitude, tmp$longitude, unit=1000) %>%
      dist() %>% max()
    tibble(max_dist=md)
  }
  
  max_dist_within_clusters = lib_meta %>%
    select(geo_cluster, latitude, longitude) %>%
    group_by(geo_cluster) %>%
    group_modify(max_dist) %>%
    ungroup()
  
  cluster_centroids = lib_meta %>%
    select(geo_cluster, latitude, longitude) %>%
    group_by(geo_cluster) %>%
    summarise(lat=mean(latitude), lon=mean(longitude), nsample=n()) %>%
    left_join(max_dist_within_clusters, by = "geo_cluster")
  
  cluster_coord = sf::st_as_sf(
    cluster_centroids,
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  pie_data = lib_meta %>%
    select(geo_cluster, genus) %>%
    group_by(geo_cluster) %>%
    summarise(
      culex=sum(genus=="Culex"), 
      anopheles=sum(genus=="Anopheles"),
      aedes=sum(genus=="Aedes"),
      armigeres=sum(genus=="Armigeres"),
      others=sum(genus%in%c("Lutzia","Mansonia","Coquillettidia","Mimomyia"))
    ) %>%
    left_join(cluster_centroids)
  
  
  ggplot() + 
    geom_sf(data=sf::st_transform(china_map, crs=4326)) + 
    scatterpie::geom_scatterpie(aes(x=lon, y=lat, r=sqrt(nsample)/10),
                                data=pie_data, cols=c("culex","anopheles","aedes","armigeres","others"), 
                                color=NA) +
    scatterpie::geom_scatterpie_legend(sqrt(pie_data$nsample)/10, x=85, y=21) +
    theme_bw() +
    ylab("Latitude") +
    xlab("Longitude") +
    ggsci::scale_fill_npg() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

draw_figure1b = function(tree_mos_cox1, lib_meta) {
  color_mapping = setNames(
    c("#3F5388", "#4DBBD6", "#E64B37", "#019F87", 
      "#F29B7F", "#F29B7F", "#F29B7F", "#F29B7F", "#F29B7F"),
    c("Armigeres", "Anopheles", "Culex", "Aedes", 
      "Mansonia", "Mimomyia", "Coquillettidia", "Lutzia", "Others")
  )
  tree = tree_mos_cox1

  species_df = tibble(
    label= tree$tip.label,
    lib_id = tree$tip.label) %>%
    left_join(lib_meta, by="lib_id") %>%
    mutate(genus = ifelse(genus %in% c("Lutzia","Mansonia","Coquillettidia","Mimomyia"), 
                          "Others", 
                          genus))
  
  d = full_join(tree, species_df)
  
  ggtree(d) + 
    geom_tree(show.legend = F) +
    geom_tippoint(aes(color=genus), size=0.5) +
    scale_color_manual(values=color_mapping) + 
    geom_treescale(x=0.4) +
    theme(legend.position = "none",
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
}

prepare_figure1c_data = function(lib_meta, interval=10) {
  sam = function(n){
    select(lib_meta, species) %>%
      sample_n(n) %>%
      `$`(species) %>%
      unique() %>%
      length()
  }
  
  rf = function(n, rep=100) {
    arr = plyr::aaply(rep(n, rep), 1, sam)
    c(mean(arr), sd(arr))
  }
  
  plyr::adply(seq(1, nrow(lib_meta), interval), 1, function(x){rf(x, 10)})
}

draw_figure1c = function(rf, interval = 10) {
  ggplot(aes(x=1+(as.integer(X1)-1)*interval, y=V1), data=rf) +
    geom_line() + 
    geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=.1, color="blue", alpha=0.15) +
    xlab("Num. of samples") +
    ylab("Num. of\nmosquito species") +
    theme_bw()
}

draw_figure1d = function(rpm_table, lib_meta, virus_meta) {
  core_virome = filter(virus_meta, is_core)$virus_name
  total_virome = virus_meta$virus_name
  
  total_virome_reads = rpm_table %>%
    select(-lib_id) %>%
    rowSums() %>%
    bind_cols(lib_id = rpm_table$lib_id, total_reads_all = .) %>%
    left_join(select(lib_meta, lib_id, norRNA_reads)) %>%
    mutate(total_reads_all = total_reads_all/1e6 * norRNA_reads) %>%
    select(-norRNA_reads)
  
  total_rpm_key = rpm_table %>%
    select(all_of(core_virome)) %>%
    rowSums() %>%
    bind_cols(lib_id = rpm_table$lib_id, total_rpm_key = .)
  
  key_virus_richness = rpm_table %>%
    select(lib_id, all_of(core_virome)) %>%
    pivot_longer(cols=2:ncol(.), names_to="virus_name", values_to="incidence") %>%
    mutate(incidence = ifelse(incidence > 0, 1, 0)) %>%
    pivot_wider(names_from="virus_name", values_from="incidence") %>%
    select(-lib_id) %>%
    rowSums() %>%
    bind_cols(lib_id = rpm_table$lib_id, key_virus_richness = .)
  
  rpm_stat = tibble(lib_id = lib_meta$lib_id) %>%
    left_join(total_virome_reads, by="lib_id") %>%
    left_join(total_rpm_key, by="lib_id") %>%
    left_join(key_virus_richness, by="lib_id") %>%
    replace_na(list(total_virome_reads = 0, 
                    total_rpm_key = 0, 
                    key_virus_richness = 0)) %>%
    left_join(select(lib_meta, lib_id, genus_alt, species, norRNA_reads), 
              by="lib_id") %>%
    arrange(genus_alt, desc(total_reads_all))  %>%
    mutate(total_reads_key = total_rpm_key /1e6 * norRNA_reads)
  
  ggplot(data=rpm_stat) +
    geom_bar(aes(x=1:nrow(rpm_stat),y=log10(norRNA_reads*2+1),
                 fill = "total non-rRNA reads"), 
             stat="identity") + 
    geom_bar(aes(x=1:nrow(rpm_stat),y=log10(total_reads_all+1),
                 fill="total-virome reads"), 
             stat="identity") +
    geom_bar(aes(x=1:nrow(rpm_stat),y=log10(total_reads_key+1),
                 fill = "core-virome reads"), 
             stat="identity") + 
    scale_fill_manual(values=c("total non-rRNA reads"="#D9D9D1", 
                        "total-virome reads"="#FCE2C4", 
                        "core-virome reads"="#92BDC8")) +
    lemon::facet_rep_grid(~genus_alt,scales="free_x", space = "free" ) + 
    ylab("log(reads+1)") +
    xlab("") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ylim(c(0,8)) + 
    lemon::coord_capped_cart(left='both') +
    
    theme(panel.background = element_rect(fill = "white"),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          panel.grid.major.y = element_line(color ="white", linewidth=0.2),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color ="black", linewidth=0.4),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text = element_text(size=10, color="black"),
          axis.title.y = element_text(size = 12))
}

draw_figure1e = function(rpm_table, lib_meta, virus_meta) {
  lib_meta = lib_meta %>%
    left_join(count(., species, name="sp_count")) %>%
    arrange(genus, desc(sp_count), species) %>%
    mutate(display_order = 1:nrow(.))
  
  virus_meta = arrange(virus_meta, display_order)
  
  core_virome = filter(virus_meta, is_core)$virus_name
  total_virome = virus_meta$virus_name
  
  virus_meta = filter(virus_meta, is_core)
  rpm_key = rpm_table %>%
    select(lib_id, all_of(core_virome)) %>%
    as.data.frame()

  rownames(rpm_key) = rpm_key$lib_id
  rpm_key = rpm_key[lib_meta$lib_id,]
  
  incidence_key_by_genus = rpm_key %>%
    left_join(select(lib_meta, lib_id, genus=genus_alt), by="lib_id") %>%
    select(-lib_id) %>%
    group_by(genus) %>%
    summarise_all(~sum(.x>0))
  
  incidence_key_by_genus_grp = incidence_key_by_genus %>%
    pivot_longer(cols=2:ncol(.), names_to="virus_name", values_to="incidence") %>%
    left_join(virus_meta) %>%
    group_by(genus, virus_group) %>%
    summarise(incidence=sum(incidence))
  
  rpm_key_by_grp = rpm_key %>%
    pivot_longer(cols=2:ncol(.), names_to="virus_name", values_to="RPM") %>%
    #filter(RPM > 0) %>%
    left_join(virus_meta) %>%
    group_by(lib_id, virus_group) %>%
    summarise(RPM=sum(RPM)) %>%
    filter(!is.na(lib_id)) %>%
    pivot_wider(names_from = "virus_group", values_from="RPM", values_fill=0) %>%
    as.data.frame()
  
  rownames(rpm_key_by_grp) = rpm_key_by_grp$lib_id
  rpm_key_by_grp$lib_id = NULL
  
  display_order = virus_meta %>% 
    group_by(baltimore_class, virus_group) %>% 
    summarise(display_order = min(display_order)) %>% 
    ungroup() %>%
    arrange(display_order) 
  
  rpm_key_by_grp = rpm_key_by_grp[,display_order$virus_group]
  lg_rpm_key_by_grp = t(log10(as.matrix(rpm_key_by_grp) + 1))
  RPM_col_fun = colorRamp2(c(0, ceiling(max(lg_rpm_key_by_grp))), c("#F4F2ED", "#DD1332"))
  
  baltimore_class = unique(display_order$baltimore_class)
  baltimore_class_color = c(brewer.pal(length(baltimore_class), "Set3"))
  row_ano = rowAnnotation(virus_taxa = anno_simple(
    display_order$baltimore_class,
    col = setNames(baltimore_class_color, baltimore_class),
    width = unit(0.3, "cm")
  ))
  
  ht = Heatmap(
    lg_rpm_key_by_grp,
    cluster_columns = F,
    cluster_rows = F,
    col = RPM_col_fun,
    row_order = 1:nrow(lg_rpm_key_by_grp),
    column_order = 1:ncol(lg_rpm_key_by_grp),
    column_split = factor(lib_meta$genus_alt, levels = c("Aedes", "Anopheles","Armigeres","Culex",  "Others"), ordered = T), 
    column_labels = rep("", ncol(lg_rpm_key_by_grp)),
    #column_title_rot = 30,
    left_annotation = row_ano,
    use_raster = T,
    name="log(RPM+1)"
  )
  lgd = Legend(
    labels = baltimore_class,
    title = "Baltimore class",
    legend_gp = gpar(fill = baltimore_class_color)
  )
  draw(ht, annotation_legend_list = list(lgd))

}

draw_figure1f = function(rpm_table, lib_meta, virus_meta) {
  lib_meta = lib_meta %>%
    left_join(count(., species, name="sp_count")) %>%
    arrange(genus, desc(sp_count), species) %>%
    mutate(display_order = 1:nrow(.))
  
  virus_meta = arrange(virus_meta, display_order)
  
  core_virome = filter(virus_meta, is_core)$virus_name
  total_virome = virus_meta$virus_name
  
  virus_meta = filter(virus_meta, is_core)
  rpm_key = rpm_table %>%
    select(lib_id, all_of(core_virome)) %>%
    as.data.frame()
  
  rownames(rpm_key) = rpm_key$lib_id
  rpm_key = rpm_key[lib_meta$lib_id,]

  incidence_key_by_genus = rpm_key %>%
    left_join(select(lib_meta, lib_id, genus=genus_alt), by="lib_id") %>%
    select(-lib_id) %>%
    group_by(genus) %>%
    summarise_all(~sum(.x>0))
  
  incidence_key_by_genus_grp = incidence_key_by_genus %>%
    pivot_longer(cols=2:ncol(.), names_to="virus_name", values_to="incidence") %>%
    left_join(virus_meta) %>%
    group_by(genus, virus_group) %>%
    summarise(incidence=sum(incidence))
  
  ggplot(data=incidence_key_by_genus_grp) + 
    geom_bar(aes(x=genus, y=incidence, fill=virus_group), stat="identity", position="fill") + 
    theme_bw() +
    theme(axis.title = element_text(size=6),
          axis.text = element_text(size=5),
          legend.text = element_text(size=5),
          legend.title = element_text(size=6),
          legend.key.size = unit(2, "mm"),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_fill_manual(values=c("#40655e", "#8ce3e6", "#3e3c8d", "#b9b6cb",
                               "#922eb1", "#f764de", "#daa4f9", "#3693f2",
                               "#5bef8f", "#1e965f", "#bde267", "#673d17",
                               "#f6a679", "#fd5917", "#9a2a06", "#e5196a",
                               "#902d54", "#ebc30e"))
  
}

draw_virus_prevalence = function(rpm_table, lib_meta, virus_meta) {
  lib_meta = lib_meta %>%
    left_join(count(., species, name="sp_count")) %>%
    arrange(genus, desc(sp_count), species) %>%
    mutate(display_order = 1:nrow(.))
  
  virus_meta = arrange(virus_meta, display_order)
  
  core_virome = filter(virus_meta, is_core)$virus_name
  total_virome = virus_meta$virus_name
  
  virus_meta = filter(virus_meta, is_core)
  rpm_key = rpm_table %>%
    select(lib_id, all_of(core_virome)) %>%
    as.data.frame()
  
  rownames(rpm_key) = rpm_key$lib_id
  rpm_key = rpm_key[lib_meta$lib_id,]
  
  rpm_key_by_grp = rpm_key %>%
    pivot_longer(cols=2:ncol(.), names_to="virus_name", values_to="RPM") %>%
    left_join(virus_meta, by="virus_name") %>%
    left_join(lib_meta, by="lib_id") %>%
    
    group_by(lib_id, virus_group) %>%
    summarise(incidence=(sum(RPM)>0)*1, species=first(species)) %>%
    group_by(species, virus_group) %>%
    summarise(prevalence=sum(incidence)/length(unique(lib_id))) %>%
    pivot_wider(names_from = "virus_group", values_from="prevalence", values_fill=0) %>%
    as.data.frame()
  
  virus_display_order = virus_meta %>% 
    group_by(baltimore_class, virus_group) %>% 
    summarise(display_order = min(display_order)) %>% 
    ungroup() %>%
    arrange(display_order)
  
  mosquito_display_order = lib_meta %>% 
    group_by(genus_alt, species) %>% 
    summarise(display_order = min(display_order)) %>% 
    ungroup() %>%
    arrange(display_order)
  
  rownames(rpm_key_by_grp) = rpm_key_by_grp$species
  rpm_key_by_grp$lib_id = NULL
  
  rpm_key_by_grp = t(rpm_key_by_grp[mosquito_display_order$species,virus_display_order$virus_group])
  col_fun = colorRamp2(c(0, ceiling(max(rpm_key_by_grp))), c("#F4F2ED", "#DD1332"))
  
  assertthat::assert_that(max(rpm_key_by_grp)==1)
  
  baltimore_class = unique(virus_display_order$baltimore_class)
  baltimore_class_color = c(brewer.pal(length(baltimore_class), "Set3"))
  row_ano = rowAnnotation(virus_taxa = anno_simple(
    virus_display_order$baltimore_class,
    col = setNames(baltimore_class_color, baltimore_class),
    width = unit(0.3, "cm")
  ))
  
  spcount = count(lib_meta, species) %>% as.data.frame()
  rownames(spcount) = spcount$species
  spcount = spcount[mosquito_display_order$species,]
  col_ano = columnAnnotation(
    num_individuals = anno_barplot(spcount$n)
  )

  colnames(rpm_key_by_grp) = str_replace(colnames(rpm_key_by_grp), "^Culex\\.? ?", "Cu. ")
  colnames(rpm_key_by_grp) = str_replace(colnames(rpm_key_by_grp), "^Aedes\\.? ?", "Ae. ")
  colnames(rpm_key_by_grp) = str_replace(colnames(rpm_key_by_grp), "^Anopheles\\.? ?", "An. ")
  colnames(rpm_key_by_grp) = str_replace(colnames(rpm_key_by_grp), "^Armigeres\\.? ?", "Ar. ")
  colnames(rpm_key_by_grp) = str_replace(colnames(rpm_key_by_grp), "^Lutzia\\.? ?", "Lu. ")
  colnames(rpm_key_by_grp) = str_replace(colnames(rpm_key_by_grp), "^Mimomyia\\.? ?", "Mi. ")
  colnames(rpm_key_by_grp) = str_replace(colnames(rpm_key_by_grp), "^Mansonia\\.? ?", "Ma. ")
  colnames(rpm_key_by_grp) = str_replace(colnames(rpm_key_by_grp), "^Coquillettidia\\.? ?", "Co. ")
  
  ht = Heatmap(
    rpm_key_by_grp,
    cluster_columns = F,
    cluster_rows = F,
    col = col_fun,
    row_order = 1:nrow(rpm_key_by_grp),
    column_order = 1:ncol(rpm_key_by_grp),
    column_split = factor(mosquito_display_order$genus_alt, levels = c("Aedes", "Anopheles","Armigeres","Culex",  "Others"), ordered = T), 
    #column_labels = rep("", ncol(rpm_key_by_grp)),
    #column_title_rot = 30,
    top_annotation = col_ano,
    left_annotation = row_ano,
    use_raster = F,
    name="Prevalence"
  )
  lgd = Legend(
    labels = baltimore_class,
    title = "Baltimore class",
    legend_gp = gpar(fill = baltimore_class_color)
  )
  draw(ht, annotation_legend_list = list(lgd))
}
