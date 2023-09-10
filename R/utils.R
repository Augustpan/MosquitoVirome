library(tidyverse)

create_output_csv = function(data, filename) {
  write_csv(data, filename, na="")
  filename
}

create_output_tsv = function(data, filename) {
  write_tsv(data, filename, na="")
  filename
}

create_output_figure = function(plot, filename, width, height, units) {
  ggsave(plot = plot, 
         filename = filename, 
         width = width, 
         height = height,
         units = units)
  filename
}

create_output_figure_CairoPDF = function(p, filename, width, height) {
  Cairo::CairoPDF(file = filename,       
                  width = width,
                  height = height)
  plot(p)
  dev.off()
  filename
}

rename_contig = function(tbl_viral_genome_cov_not_masked,
                         tbl_mask_stat) {
  name_map = setNames(tbl_mask_stat$virus_name, tbl_mask_stat$contig_name)
  tbl_viral_genome_cov_not_masked %>%
    mutate(virus_name = name_map[virus_name])
}

pad_rpm_table = function(rpm_table, lib_meta) {
  bind_rows(
    rpm_table,
    map_df(setdiff(lib_meta$lib_id, rpm_table$lib_id), 
           function(x) { setNames(c(x, rep(0, ncol(rpm_table)-1)),
                                  colnames(rpm_table))}) %>%
      mutate(across(2:ncol(.), as.double))
  )
}
