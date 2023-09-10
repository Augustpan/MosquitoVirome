library(tidyverse)

rpm_table = read_csv("output/rpm_table_virus_masked.csv")
lib_meta = read_csv("raw_data/lib_meta.csv")
vm = read_csv("output/virus_meta_merged.csv")
core_virome = filter(vm, is_core)$virus_name
rpm_table = select(rpm_table, lib_id, all_of(core_virome))

virus_names = colnames(rpm_table)[2:ncol(rpm_table)]
host_df = select(lib_meta, lib_id, genus, species)

max_rpm_per_species = rpm_table %>%
    left_join(host_df, by="lib_id") %>%
    relocate("species", .after=1) %>% 
    relocate("genus", .after=1) %>%
    select(-lib_id, -genus) %>%
    group_by(species) %>%
    summarise_all(max)

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

cst_virus = filter(virus_host_stat, num_host_species>1)$virus_name
max_rpm_cst = max_rpm_per_species %>%
    select(species, all_of(cst_virus)) %>%
    pivot_longer(cols=2:ncol(.), values_to="RPM") %>%
    filter(RPM>0)

for (virus in cst_virus) {
    nsp = length(filter(max_rpm_cst, name == "narna0001")$species)
    g= ggplot(aes(x=species, y=log10(RPM+1)), data=filter(max_rpm_cst, name == virus)) + 
        geom_bar(stat="identity") + 
        geom_hline(yintercept=3, size=0.2, linetype="dashed") + 
        geom_hline(yintercept=1, size=0.2, linetype="dashed") +
        ggtitle(virus) + 
        ylab("Max logRPM") +
        theme(axis.text.x = element_text(angle=30, hjust=1))
    ggsave(str_glue("../output/CST_RPM_stat/{virus}.jpg"), plot=g, width=0.1+0.2*nsp, height=3)
}

v = filter(max_rpm_cst, RPM>1) %>% 
    group_by(name) %>% 
    summarise(nsp=sum(RPM>0)) %>% 
    filter(nsp>1) %>%
    `$`(name)
filter(max_rpm_cst, RPM>1) %>%
    select(species, virus_name=name) %>%
    left_join(distinct(host_df, genus, species)) %>%
    relocate(genus, .before=1) %>%
    filter(virus_name %in% v) %>%
    write_csv("output/cst_all.csv")

v = filter(max_rpm_cst, RPM>10) %>% 
    group_by(name) %>% 
    summarise(nsp=sum(RPM>0)) %>% 
    filter(nsp>1) %>%
    `$`(name)
filter(max_rpm_cst, RPM>10) %>%
    select(species, virus_name=name) %>%
    left_join(distinct(host_df, genus, species)) %>%
    relocate(genus, .before=1) %>%
    filter(virus_name %in% v) %>%
    write_csv("output/cst_hiconf.csv")

v = filter(max_rpm_cst, RPM>100) %>% 
    group_by(name) %>% 
    summarise(nsp=sum(RPM>0)) %>% 
    filter(nsp>1) %>%
    `$`(name)
filter(max_rpm_cst, RPM>100) %>%
    select(species, virus_name=name) %>%
    left_join(distinct(host_df, genus, species)) %>%
    relocate(genus, .before=1) %>%
    filter(virus_name %in% v) %>%
    write_csv("output/cst_vhconf.csv")

v = filter(max_rpm_cst, RPM>1000) %>% 
    group_by(name) %>% 
    summarise(nsp=sum(RPM>0)) %>% 
    filter(nsp>1) %>%
    `$`(name)
filter(max_rpm_cst, RPM>1000) %>%
    select(species, virus_name=name) %>%
    left_join(distinct(host_df, genus, species)) %>%
    relocate(genus, .before=1) %>%
    filter(virus_name %in% v) %>%
    write_csv("output/cst_uhconf.csv")
