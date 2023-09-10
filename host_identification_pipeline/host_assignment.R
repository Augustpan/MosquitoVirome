library(tidyverse)
THRESHOLD_PIDENT = 40
THRESHOLD_LENGTH = 100

VMR = readxl::read_excel("VMR_MSL38_v1.xlsx") %>%
  separate_rows(`Virus GENBANK accession`, sep = ";\\s*") %>%
  mutate(`Virus GENBANK accession` = str_replace(`Virus GENBANK accession`, "(^.+?:\\s*)|(\\s*\\(.+\\)\\s*)", "")) %>%
  separate_rows(`Virus REFSEQ accession`, sep = ";\\s*") %>%
  mutate(`Virus REFSEQ accession` = str_replace(`Virus REFSEQ accession`, "(^.+?:\\s*)|(\\s*\\(.+\\)\\s*)", "")) %>%
  pivot_longer(cols = c(`Virus GENBANK accession`, `Virus REFSEQ accession`),
               names_to = "accession_type",
               values_to = "accession")

VHDB = read_tsv("virushostdb.tsv") %>%
  separate_rows(`refseq id`, sep=",\\s*") %>% 
  group_by(`refseq id`) %>% 
  summarise(
    across(1:(ncol(.)-1), first),
    host_tax_ids = paste(unique(`host tax id`), collapse = "/"),
    host_lineages = paste(unique(`host lineage`), collapse = "/"))

blast_table_header = c(
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
)

blast_vs_vmr = read_tsv("query_vs_vmr.txt", col_names = blast_table_header)

sseqid_parsed = str_extract(blast_vs_vmr$sseqid,
            "lcl\\|(.+).\\d+_prot_(.+).\\d+_\\d+",
            group = c(1, 2))
blast_vs_vmr$sseqid_nucl_acc = sseqid_parsed[,1]
blast_vs_vmr$sseqid_prot_acc = sseqid_parsed[,2]
blast_vs_vmr = left_join(blast_vs_vmr, 
                          select(VMR, accession, `Host source`), 
                          by = c("sseqid_nucl_acc"="accession")) %>%
  mutate(vmr_invert_assoc = str_detect(`Host source`, "invertebrates"))

blast_vs_vhdb = read_tsv("query_vs_vhdb.txt", col_names = blast_table_header)
sseqid_parsed = str_extract(blast_vs_vhdb$sseqid,
            "^(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*)",
            group = c(1, 4))
blast_vs_vhdb$sseqid_prot_acc = sseqid_parsed[,1]
blast_vs_vhdb$sseqid_nucl_acc = sseqid_parsed[,2]
blast_vs_vhdb = left_join(blast_vs_vhdb, 
          select(VHDB, `refseq id`, host_lineages), 
          by = c("sseqid_nucl_acc"="refseq id")) %>%
  mutate(vhdb_invert_assoc = ifelse(host_lineages != "NA", str_detect(host_lineages, "Arthropoda"), NA))

blast_vs_all = read_tsv("all_vs_all.txt", col_names = blast_table_header)

manual_assignment = read_csv("manual_assignment.csv") %>%
  select(qseqid, assignment_rpm, assignment_manual) %>%
  mutate(assignment_rpm = ifelse(assignment_rpm==F, NA, T))

vhdb_info = blast_vs_vhdb %>%
  filter(pident > THRESHOLD_PIDENT, 
         length > THRESHOLD_LENGTH) %>%
  arrange(bitscore) %>%
  group_by(qseqid) %>%
  summarise(vhdb_invert_assoc_score = sum(vhdb_invert_assoc),
            vhdb_invert_assoc_tophit = first(vhdb_invert_assoc),
            vhdb_tophit_pident = first(pident),
            vhdb_tophit_evalue = first(evalue),
            vhdb_tophit_length = first(length))

vmr_info = blast_vs_vmr %>%
  filter(pident > THRESHOLD_PIDENT, 
         length > THRESHOLD_LENGTH) %>%
  arrange(bitscore) %>%
  group_by(qseqid) %>%
  summarise(vmr_invert_assoc_score = sum(vmr_invert_assoc),
            vmr_invert_assoc_tophit = first(vmr_invert_assoc),
            vmr_tophit_pident = first(pident),
            vmr_tophit_evalue = first(evalue),
            vmr_tophit_length = first(length))

assign_host = function(assignment_vmr, 
                       assignment_vhdb, 
                       assignment_manual, 
                       assignment_rpm) {
  
  if (!is.na(assignment_manual))
    return(paste0(assignment_manual, "_MANUAL"))
  
  if (!is.na(assignment_vmr))
    return(paste0(assignment_vmr, "_VMR"))
  
  if (!is.na(assignment_rpm))
    return(paste0(assignment_rpm, "_RPM"))
  
  if (!is.na(assignment_vhdb))
    return(paste0(assignment_vhdb, "_VHDB"))

  return(NA)
}

host_table = tibble(qseqid = unique(blast_vs_all$qseqid)) %>%
  left_join(vmr_info, by = "qseqid") %>%
  left_join(vhdb_info, by = "qseqid") %>%
  left_join(manual_assignment, by = "qseqid") %>%
  mutate(
    assignment_vmr = vmr_invert_assoc_tophit,
    assignment_vhdb = vhdb_invert_assoc_tophit) %>%
  rowwise() %>%
  mutate(
    host_assignment_evidence = assign_host(
      assignment_vmr, 
      assignment_vhdb, 
      assignment_manual, 
      assignment_rpm)) %>%
  separate_wider_delim(host_assignment_evidence,"_", names=c("is_core", "evidence")) %>%
  mutate(is_core = as.logical(is_core))

iteration = 1
while (TRUE) {
  propagate_table = blast_vs_all %>%
    filter(pident > THRESHOLD_PIDENT,
           length > THRESHOLD_LENGTH) %>%
    select(qseqid, sseqid) %>%
    left_join(select(host_table, sseqid=qseqid, is_core), by="sseqid") %>%
    filter(!is.na(is_core)) %>%
    group_by(qseqid) %>%
    summarise(core_hits = sum(is_core)/n()) %>%
    mutate(pc = core_hits > 0.8, 
           pnc = core_hits < 0.2,
           assignment_propagated = ifelse(pc, TRUE, ifelse(pnc, FALSE, NA)),
           evidence_propagated = ifelse(pc|pnc, str_glue("PROPAGATED_ITER{iteration}"), NA)) %>%
    select(-core_hits, -pc, -pnc) %>%
    filter(!is.na(assignment_propagated)) %>%
    filter(qseqid %in% filter(host_table, is.na(is_core))$qseqid)
  
  if (nrow(propagate_table) == 0) {
    print(str_glue("stop at iteration {iteration-1}"))
    break
  }
  
  host_table = host_table %>%
    left_join(propagate_table, by="qseqid") %>%
    mutate(to_assign = is.na(is_core),
           is_core = ifelse(to_assign, assignment_propagated, is_core),
           evidence = ifelse(to_assign, evidence_propagated, evidence)) %>%
    select(-to_assign, -assignment_propagated, -evidence_propagated)
  
  print(str_glue("{nrow(propagate_table)} records propagated"))
  iteration = iteration + 1
}

manifest = read_tsv("virus_manifest.txt", col_names=c("qseqid"))
write_excel_csv(host_table, "OUTPUT_HOST_TABLE_WITH_REF.csv")
write_excel_csv(filter(host_table, qseqid %in% manifest$qseqid), 
                "OUTPUT_HOST_TABLE.csv")
