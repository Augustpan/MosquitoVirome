# _targets.R file

library(targets)
library(tidyverse)

source("R/utils.R")
source("R/rpm_table_cleanning.R")
source("R/Figure1.R")
source("R/Figure2.R")

tar_option_set(packages = c("tidyverse", "lemon", "sf", "ggtree",
                            "RColorBrewer", "circlize", "ComplexHeatmap"))

list(
  #### generating RPM table ####
  
  # define raw data filenames
  tar_target(file_viral_genome_cov_all_masked, 
             "raw_data/viral_genomes_mapping/cov_all_rep_seqs_all_masked.txt", 
             format = "file"),
  tar_target(file_viral_genome_cov_host_masked, 
             "raw_data/viral_genomes_mapping/cov_all_rep_seqs_host_masked.txt", 
             format = "file"),
  tar_target(file_viral_genome_cov_virus_masked, 
             "raw_data/viral_genomes_mapping/cov_all_rep_seqs_reciprocal_masked.txt",
             format = "file"),
  tar_target(file_viral_genome_cov_not_masked, 
             "raw_data/viral_genomes_mapping/cov_all_rep_seqs.txt", 
             format = "file"),
  tar_target(file_lane_info, 
             "raw_data/lane_info.txt", 
             format = "file"),
  tar_target(file_norrna_reads, 
             "raw_data/norRNA_reads.txt", 
             format = "file"),
  tar_target(file_mask_stat, 
             "raw_data/viral_genomes/mask_stat.txt", 
             format = "file"),
  
  # load raw data files
  tar_target(tbl_mask_stat,
             read_tsv(file_mask_stat)),
  tar_target(tbl_viral_genome_cov_all_masked, 
             load_coverage_file(file_viral_genome_cov_all_masked)),
  tar_target(tbl_viral_genome_cov_host_masked_raw, 
             load_coverage_file(file_viral_genome_cov_host_masked)),
  tar_target(tbl_viral_genome_cov_virus_masked_raw, 
             load_coverage_file(file_viral_genome_cov_virus_masked)),
  tar_target(tbl_viral_genome_cov_not_masked_raw, 
             load_coverage_file(file_viral_genome_cov_not_masked)),
  
  # rename representative contig to virus_name
  tar_target(tbl_viral_genome_cov_virus_masked,
             rename_contig(tbl_viral_genome_cov_virus_masked_raw, 
                           tbl_mask_stat)),
  tar_target(tbl_viral_genome_cov_host_masked,
             rename_contig(tbl_viral_genome_cov_host_masked_raw, 
                           tbl_mask_stat)),
  tar_target(tbl_viral_genome_cov_not_masked,
             rename_contig(tbl_viral_genome_cov_not_masked_raw, 
                           tbl_mask_stat)),
  tar_target(tbl_lane_info, 
             read_tsv(file_lane_info)),
  tar_target(tbl_norrna_reads, 
             mutate(read_tsv(file_norrna_reads), 
                    norRNA_reads=norRNA_reads*2)),
  
  # calculate RPM tables
  tar_target(rpm_table_all_masked_raw,
             calc_rpm(tbl_viral_genome_cov_all_masked, 
                      tbl_lane_info,
                      tbl_norrna_reads)),
  tar_target(rpm_table_all_masked,
             pad_rpm_table(rpm_table_all_masked_raw, 
                           tbl_lib_meta)),
  
  tar_target(rpm_table_host_masked_raw,
             calc_rpm(tbl_viral_genome_cov_host_masked, 
                      tbl_lane_info,
                      tbl_norrna_reads)),
  tar_target(rpm_table_host_masked,
             pad_rpm_table(rpm_table_host_masked_raw, 
                           tbl_lib_meta)),
  
  tar_target(rpm_table_virus_masked_raw, 
             calc_rpm(tbl_viral_genome_cov_virus_masked, 
                      tbl_lane_info,
                      tbl_norrna_reads)),
  tar_target(rpm_table_virus_masked,
             pad_rpm_table(rpm_table_virus_masked_raw, 
                           tbl_lib_meta)),
  
  tar_target(rpm_table_not_masked_raw, 
             calc_rpm(tbl_viral_genome_cov_not_masked, 
                      tbl_lane_info,
                      tbl_norrna_reads)),
  tar_target(rpm_table_not_masked,
             pad_rpm_table(rpm_table_not_masked_raw, 
                           tbl_lib_meta)),
  
  # write RPM tables to file
  tar_target(output_rpm_table_all_masked,
             create_output_csv(rpm_table_all_masked, 
                               "output/rpm_table_all_masked.csv"),
             format = "file"),
  tar_target(output_rpm_table_host_masked,
             create_output_csv(rpm_table_host_masked,
                               "output/rpm_table_host_masked.csv"),
             format = "file"),
  tar_target(output_rpm_table_virus_masked,
             create_output_csv(rpm_table_virus_masked,
                               "output/rpm_table_virus_masked.csv"),
             format = "file"),
  tar_target(output_rpm_table_not_masked,
             create_output_csv(rpm_table_not_masked, 
                               "output/rpm_table_not_masked.csv"),
             format = "file"),
  
  #### virus_meta and lib_meta ####
  # define raw data filenames
  tar_target(file_virus_host_annot,
             "raw_data/host_taxa_annotations.csv",
             format = "file"),
  tar_target(file_virus_taxa_annot,
             "raw_data/virus_taxa_annotations.csv",
             format = "file"),
  tar_target(file_lib_meta,
             "raw_data/lib_meta.csv",
             format = "file"),
  
  # load raw data files
  tar_target(tbl_virus_host_annot,
             mutate(read_csv(file_virus_host_annot), 
                    is_core = ifelse(is.na(manual_correction), 
                                     is_core, 
                                     manual_correction))),
  tar_target(tbl_virus_taxa_annot,
             read_csv(file_virus_taxa_annot)),

  # prepare virus_meta and lib_meta
  tar_target(tbl_lib_meta,
             read_csv(file_lib_meta) %>%
               mutate(genus_alt = ifelse(genus %in% c("Anopheles","Culex","Armigeres","Aedes"),
                                         genus,
                                         "Others"),
                      geo_cluster = SoDA::geoXY(latitude, 
                                                longitude, 
                                                unit = 1000) %>%
                        dist() %>%
                        hclust(method = "average") %>%
                        cutree(h = 100) %>%
                        paste0("c", .)) %>%
               left_join(tbl_lane_info, by="lib_id") %>%
               left_join(tbl_norrna_reads, by="lib_id")),
  
  tar_target(tbl_virus_meta,
             full_join(tbl_virus_taxa_annot, 
                       tbl_virus_host_annot,
                       by = "virus_name") %>%
               left_join(tbl_mask_stat, by="virus_name")),
  
  ## write to files
  tar_target(output_virus_meta,
             create_output_csv(tbl_virus_meta, 
                               "output/virus_meta_merged.csv"),
             format = "file"),
  tar_target(output_lib_meta,
             create_output_csv(tbl_lib_meta, 
                               "output/lib_meta_merged.csv"),
             format = "file"),
  
  #### Figure 1 ####
  
  # Figure1a: map
  tar_target(file_china_map,
             "map_data/Province.shp", 
             format = "file"),
  tar_target(sf_china_map, 
             sf::read_sf(file_china_map)),
  tar_target(plot_figure_1a,
             draw_figure1a(tbl_lib_meta, 
                           sf_china_map)),
  tar_target(output_figure_1a,
             create_output_figure(plot_figure_1a,
                                  "output/Figure1a.pdf",
                                  width = 8, height = 5,
                                  units = "in"),
             format = "file"),
  
  # Figure1b: mosquito phylogenetic tree
  tar_target(file_mos_tree,
             "raw_data/non_contaminated_COI.trimal.aln.treefile", 
             format = "file"),
  tar_target(tree_mos_cox1, 
             read.tree(file_mos_tree) %>%
               phytools::midpoint.root()),
  tar_target(plot_figure_1b,
             draw_figure1b(tree_mos_cox1, tbl_lib_meta)),
  tar_target(output_figure_1b,
             create_output_figure(plot_figure_1b,
                                  "output/Figure1b.pdf",
                                  width = 35, height = 45,
                                  units = "mm"),
             format = "file"),
  
  # Figure1c: mosquito rarefaction
  tar_target(data_figure_1c,
             prepare_figure1c_data(tbl_lib_meta)),
  tar_target(plot_figure_1c,
             draw_figure1c(data_figure_1c)),
  tar_target(output_figure_1c,
             create_output_figure(plot_figure_1c,
                                  "output/Figure1c.pdf",
                                  width = 79, height = 56,
                                  units = "mm"),
             format = "file"),
  
  # Figure1d: viral RNA
  tar_target(plot_figure_1d,
             draw_figure1d(rpm_table_virus_masked, 
                          tbl_lib_meta,
                          tbl_virus_meta)),
  tar_target(output_figure_1d,
             create_output_figure(plot_figure_1d,
                                  "output/Figure1d.pdf",
                                  width = 12, height = 2,
                                  units = "in"),
             format = "file"),
  
  # Figure1e: heatmap
  tar_target(plot_figure_1e,
             draw_figure1e(rpm_table_virus_masked, 
                           tbl_lib_meta,
                           tbl_virus_meta)),
  tar_target(output_figure_1e,
             create_output_figure_CairoPDF(plot_figure_1e,
                                           "output/Figure1e.pdf",
                                           width = 12, height = 3.5),
             format = "file"),
  
  # Figure1f: prevalence
  tar_target(plot_figure_1f,
             draw_figure1f(rpm_table_virus_masked, 
                           tbl_lib_meta,
                           tbl_virus_meta)),
  tar_target(output_figure_1f,
             create_output_figure(plot_figure_1f,
                                  "output/Figure1f.pdf",
                                  width = 50, height = 33,
                                  units = "mm"),
             format = "file"),
  
  #### Figure 2 ####
  tar_target(file_refseq_annot, 
             "raw_data/refseq_host_annotation.csv", 
             format = "file"),
  tar_target(tbl_refseq_annot, 
             read_csv(file_refseq_annot)),
  tar_target(plot_figure2,
             draw_figure2(tbl_virus_meta,
                          tbl_refseq_annot,
                          rpm_table_virus_masked)),
  tar_target(output_figure2,
             create_output_figure(plot_figure2,
                                  "output/Figure2.pdf",
                                  width = 160, height = 154,
                                  units = "mm"),
             format = "file")
  
  #### Figure 3 ####
  
)
