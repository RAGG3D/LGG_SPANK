list_TCGA <- c("CHOL", "UCS", "DLBC", "UVM", "MESO", "ACC", "KICH", "TGCT", "READ", 
               "PCPG", "PAAD", "SARC", "KIRP", "LIHC", "BLCA", "STAD", "COAD", "SKCM",
               "PRAD", "LUSC", "THCA", "LGG", "HNSC", "KIRC", "UCEC", "LUAD", "OV", "GBM",
               "CESC", "ESCA", "LAML", "THYM", "BRCA")

no0_CIBERSORT <- function(cancer){
  cell <- read_csv(paste0("data/TCGA_", cancer, "_CIBERSORT_PDGFDD.csv")) %>%
    tidybulk::rename(sample = `Input Sample`) %>%
    dplyr::select(-`P-value`, -`Pearson Correlation`, -RMSE) %>%
    gather(celltype, fraction, -sample)
  no0cell <- foreach(i = as.character(as.data.frame(cell %>% group_by(sample) %>% summarise() %>% ungroup())[,1]), .combine = bind_rows) %do% {
    n = cell %>% filter(sample == i) %>% arrange(fraction)
    if(as.numeric(n[1,3]) > 0) {num = 1} else{
      foreach (j = 1:nrow(n)) %do% {
        if(as.numeric(n[j,3]) > 0 && n[j-1, 3] <= 0)
        {num = j}
      }
      n[1:num,][3] = n[num,3]
      n
    }} %>%
    mutate(scale = log(fraction)) %>% 
    dplyr::select(-fraction)
  no0cell
}


LGG_transcript <- function(){
  foreach(i = list.files("data/Transcript/"), .combine = bind_rows) %do% {
    read_table2(paste0("data/Transcript/", i), col_names = F) %>%
      mutate(sample = i)
  }  %>% 
    tidybulk::rename(ensembl = `X1`) %>%
    tidybulk::rename(raw_count = `X2`) %>%
    separate(ensembl, c("ensembl_id", "c"), sep = "\\.") %>%
    inner_join(toTable(org.Hs.egENSEMBL)) %>%
    inner_join(toTable(org.Hs.egSYMBOL)) %>%
    dplyr::select(sample, symbol, raw_count) %>%
    mutate(sample = paste0(substr(basename(as.character(sample)), start = 1, stop = 49), ".gz")) %>%
    as.tibble() %>%
    inner_join(
      read_csv(paste0("data/gdc_sample_sheet.csv"))%>% mutate(sample = `File Name`)) %>%
    mutate(sample = `Case ID`) %>%
    dplyr::select(sample, symbol, raw_count) %>%
    tidybulk::aggregate_duplicates(.sample = sample, .abundance = raw_count, .transcript = symbol) %>%
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol) 
}

clinical_combine <- function(x, cancer) {
  x %>% 
    inner_join(read.csv(paste0("data/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))}

Tcell <- c("t_CD4_memory",  "t_CD8_memory",  "t_CD8_naive",   "t_gamma_delta", "t_helper")
NK <- c("nk_resting", "nk_primed_IL2", "nk_primed_IL2_PDGFD")
Bar_TCGAcelltype <- function(cancer){
  x <- read_csv(paste0("data/TCGA_", cancer, "_CIBERSORT_PDGFDD.csv")) %>%
    mutate(sample = `Input Sample`) %>%
    dplyr::select(sample, contains("nk_"), -`Input Sample`) %>%
    gather(nkstate, fraction, -sample) 
  y <- x %>%
    group_by(sample) %>% 
    summarise(order = sum(fraction)) %>%
    right_join(x, by = c("sample")) %>%
    tidybulk::rename(profile = fraction) %>%
    mutate(cat = "Fraction") %>%
    rbind(x %>% 
            group_by(sample) %>% 
            summarise(sum = sum(fraction)) %>%
            right_join(x) %>%
            filter(sum > 0) %>%
            mutate(profile = fraction/sum) %>%
            left_join(x %>%
                        group_by(sample) %>% 
                        summarise(sum = sum(fraction)) %>%
                        right_join(x) %>%
                        filter(sum > 0) %>%
                        mutate(profile = fraction/sum) %>%
                        filter(nkstate == "nk_primed_IL2_PDGFD") %>%
                        mutate(order = profile) %>%
                        dplyr::select(sample, order)) %>%
            dplyr::select(-sum, -fraction) %>%
            mutate(cat = "Percentage"))
} 

KM_TCGAcelltype_median <- function(i){
  x <- as.data.frame(read_csv(paste0("data/TCGA_", i, "_CIBERSORT_PDGFDD.csv")) %>%
                       mutate(sample = `Input Sample`) %>%
                       dplyr::select(sample, contains("nk_"), -`Input Sample`) %>%
                       mutate(nk_primed_IL2 = ifelse(
                         nk_primed_IL2 > median(nk_primed_IL2), "H", "L"
                       )) %>%
                       mutate(nk_primed_IL2_PDGFD = ifelse(
                         nk_primed_IL2_PDGFD > median(nk_primed_IL2_PDGFD), "H", "L"
                       )) %>% 
                       mutate(nk_resting = ifelse(
                         nk_resting > median(nk_resting), "H", "L"
                       )) %>%
                       gather(celltype, fraction, -sample) %>%
                       inner_join(read.csv(paste0("data/clinical_", tolower(i), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
                       mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
                       mutate(na = is.na(total_living_days)) %>%
                       mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
                       dplyr::select(sample, celltype, fraction, vital_status, total_living_days, age) %>%
                       mutate(vital_status = ifelse(vital_status == "Dead", 1, 0)))
}

Gene_heatmap_log10 <- function(cancer, gene){
  x <- foreach(i = gene, .combine = bind_rows) %do%{
    LGG %>% 
      filter(symbol == i) %>%
      mutate(item = log(raw_count_scaled + 0.1), cat = symbol) 
  } %>%
    dplyr::select(sample, cat, item) %>%
    spread(cat, item)
  a <- melt(cor(x %>% dplyr::select(-sample))[-c(3:4), c(3:4)])
  a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
  a <- a %>% inner_join(as.data.frame(a.p$P[-c(3:4), c(3:4)]) %>% 
                          mutate(X1 = rownames(as.data.frame(a.p$P[-c(3:4), c(3:4)]))) %>% 
                          gather(X2, p.value, -X1)) 
  
  a$stars <- cut(a$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
  
  ggplot(data = a , aes(x=`X2`, y=`X1`, fill=value)) + 
    geom_tile() + 
    geom_text(aes(label=stars)) + 
    scale_fill_viridis(limits = c(-0.2, 0.5)) +
    theme_bw(base_size = 16, base_family = "serif") +
    labs(fill = "Correlation", x = "", y = "") +
    theme(text=element_text(size=16, family="sans"),
          legend.position = "bottom", legend.key.width = unit(3, "line")) 
}

Gene_correlation <- function(cancer, gene){
  x <- foreach(i = gene, .combine = bind_rows) %do%{
    LGG %>% 
      filter(symbol == i) %>%
      mutate(item = log(raw_count_scaled + 0.1), cat = symbol) 
  } %>%
    dplyr::select(sample, cat, item) %>%
    spread(cat, item)
  a <- melt(cor(x %>% dplyr::select(-sample))[-c(3:4), c(3:4)])
  a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
  a <- a %>% inner_join(as.data.frame(a.p$P[-c(3:4), c(3:4)]) %>% 
                          mutate(X1 = rownames(as.data.frame(a.p$P[-c(3:4), c(3:4)]))) %>% 
                          gather(X2, p.value, -X1)) 
  
  a$stars <- cut(a$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
  
  a
}


PDGFsurvival <- function(cancer, n) {
  x <- LGG %>% 
    filter(symbol == "PDGFD"|symbol == "PDGFRB") %>%
    dplyr::select(sample, symbol, raw_count_scaled) %>%
    spread(symbol, raw_count_scaled) %>%
    mutate(PDGFD = factor(Hmisc::cut2(PDGFD, g = n), labels = c(1:nlevels(Hmisc::cut2(PDGFD, g = n))))) %>%
    mutate(PDGFRB = factor(Hmisc::cut2(PDGFRB, g = n), labels = c(1:nlevels(Hmisc::cut2(PDGFRB, g = n))))) }

TCGA_clinical <- function(x, cancer) {
  x %>% 
    inner_join(read.csv(paste0("data/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    dplyr::select(sample, cat, item, vital_status, total_living_days, age) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))}

Gene_Cell <- function(cancer, gene, cells){
  foreach(i = cells, .combine = bind_rows) %do%{
    LGG %>%
      filter(symbol == gene) %>%
      dplyr::select(sample, symbol, raw_count_scaled) %>%
      spread(symbol, raw_count_scaled) %>%
      mutate(gene = factor(Hmisc::cut2(!!as.name(gene), g = 2), labels = c(1:2))) %>%
      inner_join(read_csv(paste0("data/TCGA_", cancer, "_CIBERSORT_PDGFDD.csv")) %>%
                   tidybulk::rename("sample" = "Input Sample") %>%
                   mutate(cell = factor(Hmisc::cut2(!!as.name(i), g = 2), labels = c(1:2))) %>%
                   dplyr::select(sample, cell)) %>%
      unite("item", gene, cell, sep = "/", remove = T) %>%
      dplyr::select(sample, item) %>% 
      mutate(cat = paste0(gene, "/", i))
  } %>% TCGA_clinical(cancer)
}


Cell_clinical <- function(cancer, cells, n) {
  foreach(i = cells, .combine = bind_rows) %do% {
    read_csv(paste0("data/TCGA_", cancer, "_CIBERSORT_PDGFDD.csv")) %>%
      tidybulk::rename("sample" = "Input Sample") %>%
      mutate(item = factor(Hmisc::cut2(!!as.name(i), g = n), labels = c(1:nlevels(Hmisc::cut2(!!as.name(i), g = n))))) %>%
      dplyr::select(sample, item) %>%
      mutate(cat = i)
  }
}


NK_T_cor <- function(cancer) {
  x <- foreach(i = 1:length(NK), .combine = bind_rows) %do% {
    foreach(j = 1:length(Tcell), .combine = bind_rows) %do% {
      read_csv(paste0("data/TCGA_", cancer, "_CIBERSORT_PDGFDD.csv")) %>%
        tidybulk::rename(sample = `Input Sample`) %>%
        mutate(nk = !!as.name(NK[i]), t = !!as.name(Tcell[j])) %>%
        mutate(nk = ifelse(nk > median(nk), "H", "L")) %>%
        mutate(t = ifelse(t > median(t), "H", "L")) %>%
        dplyr::select(sample, nk, t) %>%
        unite(nk_t, c(nk, t), remove = T, sep = "/") %>%
        mutate(NK = NK[i], Tcell = Tcell[j]) %>%
        unite(NKT, c(NK, Tcell), remove = T, sep = "/") %>%
        inner_join(read.csv(paste0("data/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
        mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
        mutate(na = is.na(total_living_days)) %>%
        mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
        dplyr::select(sample, nk_t, NKT, vital_status, total_living_days, age) %>%
        mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))
    }
  }
  
  x <- as.data.frame(x)
  x$NKT = factor(x$NKT, levels=c("nk_resting/t_helper", "nk_resting/t_CD8_naive", "nk_resting/t_gamma_delta", 'nk_resting/t_CD4_memory', "nk_resting/t_CD8_memory",
                                 "nk_primed_IL2/t_helper", "nk_primed_IL2/t_CD8_naive", "nk_primed_IL2/t_gamma_delta", 'nk_primed_IL2/t_CD4_memory', "nk_primed_IL2/t_CD8_memory",
                                 "nk_primed_IL2_PDGFD/t_helper", "nk_primed_IL2_PDGFD/t_CD8_naive", "nk_primed_IL2_PDGFD/t_gamma_delta", 'nk_primed_IL2_PDGFD/t_CD4_memory', "nk_primed_IL2_PDGFD/t_CD8_memory"), ordered=TRUE)
  x$nk_t <- factor(x$nk_t, levels=c("L/L", "L/H", "H/L", "H/H"), ordered=TRUE)
  
  ggsurvplot(
    survfit(
      Surv(total_living_days, vital_status) ~ nk_t,
      data = x 
    ),
    data = x,
    facet.by = c("NKT"),
    conf.int = T,
    risk.table = F,
    pval = F,
    nrow = 3,
    conf.int.alpha = 0.15,
    legend.title = "NK/T Cell Fraction",
    short.panel.labs = T,
    panel.labs = list(NKT = c("ReNK/Helper T", "ReNK/Naive CD8 T", "ReNK/?æ? T", "ReNK/Memory CD4 T", "ReNK/Memory CD8 T",
                              "IL2NK/Helper T", "IL2NK/Naive CD8 T", "IL2NK/?æ? T", "IL2NK/Memory CD4 T", "IL2NK/Memory CD8 T",
                              "SPANK/Helper T", "SPANK/Naive CD8 T", "SPANK/?æ? T", "SPANK/Memory CD4 T", "SPANK/Memory CD8 T")),
    palette = "jco"
  ) + theme_bw(base_size = 15) + 
    guides(linetype = FALSE) 
  
}

SPANK_T_cor <- function(cancer) {
  x <- foreach(i = 3, .combine = bind_rows) %do% {
    foreach(j = 1:length(Tcell), .combine = bind_rows) %do% {
      read_csv(paste0("data/TCGA_", cancer, "_CIBERSORT_PDGFDD.csv")) %>%
        tidybulk::rename(sample = `Input Sample`) %>%
        mutate(nk = !!as.name(NK[i]), t = !!as.name(Tcell[j])) %>%
        mutate(nk = ifelse(nk > median(nk), "H", "L")) %>%
        mutate(t = ifelse(t > median(t), "H", "L")) %>%
        dplyr::select(sample, nk, t) %>%
        unite(nk_t, c(nk, t), remove = T, sep = "/") %>%
        mutate(NK = NK[i], Tcell = Tcell[j]) %>%
        unite(NKT, c(NK, Tcell), remove = T, sep = "/") %>%
        inner_join(read.csv(paste0("data/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
        mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
        mutate(na = is.na(total_living_days)) %>%
        mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
        dplyr::select(sample, nk_t, NKT, vital_status, total_living_days, age) %>%
        mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))
    }
  } 
  
  x <- as.data.frame(x)
  x$NKT = factor(x$NKT, levels=c("nk_primed_IL2_PDGFD/t_helper", "nk_primed_IL2_PDGFD/t_CD8_naive", "nk_primed_IL2_PDGFD/t_gamma_delta", 'nk_primed_IL2_PDGFD/t_CD4_memory', "nk_primed_IL2_PDGFD/t_CD8_memory"), ordered=TRUE)
  x$nk_t <- factor(x$nk_t, levels=c("L/L", "L/H", "H/L", "H/H"), ordered=TRUE)
  
  ggsurvplot(
    survfit(
      Surv(total_living_days, vital_status) ~ nk_t,
      data = x 
    ),
    data = x,
    facet.by = c("NKT"),
    conf.int = T,
    risk.table = F,
    pval = F,
    nrow = 1,
    conf.int.alpha = 0.15,
    legend.title = "NK/T Cell Fraction",
    short.panel.labs = T,
    panel.labs = list(NKT = c("SPANK/Helper T", "SPANK/Naive CD8 T", "SPANK/?æ? T", "SPANK/Memory CD4 T", "SPANK/Memory CD8 T")),
    palette = "jco"
  ) + theme_bw(base_size = 15) + 
    guides(linetype = FALSE) 
  
}

genelist_cancer <- function(y, cancer, genelist) {
  x <- foreach(i = genelist, .combine = bind_rows) %do% {
    y %>% filter(symbol == i) %>%
      mutate(median = factor(Hmisc::cut2(raw_count_scaled, g = 2), labels = c(1:2))) %>%
      mutate(quantile = factor(Hmisc::cut2(raw_count_scaled, g = 4), labels = c(1:nlevels(Hmisc::cut2(raw_count_scaled, g = 4))))) 
  } %>% 
    dplyr::select(sample, symbol, median, quantile) %>%
    gather(cat, item, -c(sample, symbol)) %>%
    inner_join(read.csv(paste0("data/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0)) %>%
    dplyr::select(sample, symbol, cat, item, vital_status, total_living_days, age)
  
  x <- as.data.frame(x)
  x$item = factor(x$item, levels=c(1,2,3,4))
  x$cat <- factor(x$cat, levels=c("median", "quantile"), ordered=TRUE)
  x
}

NK_T_receptor_heatmap_log10 <- function(cancer){
  Immunereceptor <- list("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "CRTAM", 
                         "CD226", "CD244", "KIR2DL4", "NCR1", "NCR2", "NCR3")
  x <- LGGno0 %>%
    spread(celltype, scale) %>%
    dplyr::select(sample, contains("nk_"), contains("t_"), -mast_cell) %>%
    inner_join(foreach(i = Immunereceptor, .combine = bind_rows) %do%{
      LGG %>% 
        filter(symbol == i) %>%
        mutate(item = log(raw_count_scaled + 1), cat = symbol) 
    } %>%
      dplyr::select(sample, cat, item) %>%
      spread(cat, item))
  
  a <- melt(cor(x %>% dplyr::select(-sample))[-c(1:8), c(1:8)])
  a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
  a <- a %>% inner_join(as.data.frame(a.p$P[-c(1:8), c(1:8)]) %>% 
                          mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:8), c(1:8)]))) %>% 
                          gather(X2, p.value, -X1)) 
  
  a$stars <- cut(a$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
  a$X2 <- factor(a$X2, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD', "t_helper", "t_CD8_naive", "t_gamma_delta", 't_CD4_memory', "t_CD8_memory"), ordered=TRUE)
  
  ggplot(data = a , aes(x=`X2`, y=`X1`, fill=value)) + 
    geom_tile() + 
    geom_text(aes(label=stars)) + 
    scale_fill_viridis() +
    labs(y = "", x = "") + #y = "Log Transformed Normalized Count", x = "Log Transformed NK Proportion"
    theme_bw() +
    scale_x_discrete(labels = c("ReNK", "IL2NK", "SPANK", "Helper\nT", "Naive\nCD8 T", "?æ? T", "Memory\nCD4 T", "Memory\nCD8 T")) +
    theme(text=element_text(size=16, family="sans"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = -30),
          legend.key.width = unit(3, "line")) +
    labs(fill = "Correlation")
}

NK_selectT_receptor_heatmap_log10 <- function(cancer){
  cancer = "LGG"
  Immunereceptor <- list("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "CRTAM", 
                         "CD226", "CD244", "KIR2DL4", "NCR1", "NCR2", "NCR3")
  x <- LGGno0 %>%
    spread(celltype, scale) %>%
    dplyr::select(sample, contains("nk_"), contains("t_"), -mast_cell) %>%
    inner_join(foreach(i = Immunereceptor, .combine = bind_rows) %do%{
      LGG %>% 
        filter(symbol == i) %>%
        mutate(item = log(raw_count_scaled + 1), cat = symbol) 
    } %>%
      dplyr::select(sample, cat, item) %>%
      spread(cat, item)) %>%
    dplyr::select(-t_gamma_delta, -t_CD8_naive, -t_CD4_memory)
  
  a <- melt(cor(x %>% dplyr::select(-sample))[-c(1:5), c(1:5)])
  a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
  a <- a %>% inner_join(as.data.frame(a.p$P[-c(1:5), c(1:5)]) %>% 
                          mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:5), c(1:5)]))) %>% 
                          gather(X2, p.value, -X1)) 
  
  a$stars <- cut(a$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
  a$X2 <- factor(a$X2, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD', "t_helper", "t_CD8_memory"), ordered=TRUE)
  
  ggplot(data = a , aes(x=`X2`, y=`X1`, fill=value)) + 
    geom_tile() + 
    geom_text(aes(label=stars)) + 
    scale_fill_viridis() +
    labs(y = "", x = "") + #y = "Log Transformed Normalized Count", x = "Log Transformed NK Proportion"
    theme_bw() +
    scale_x_discrete(labels = c("ReNK", "IL2NK", "SPANK", "Helper\nT", "Memory\nCD8 T")) +
    theme(text=element_text(size=16, family="sans"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = -30),
          legend.key.width = unit(3, "line")) +
    labs(fill = "Correlation")
}

NK_T_genelist_heatmap_log10 <- function(cancer, genelist){
  cancer = "LGG"
  a = 42
  b = 51
  c = 54
  d = 61
  genelist <- list("CCNE2", "CDC20", "ZWINT", "UHRF1", "ERCC6L", "CDCA3", "PLK1", "PRC1", "CCNF", "CCNA2",
                   "CENPE", "KIFC1", "MAD2L1", "CDK1", "CCNB1", "NUP37", "CEP55", "CHTF18", "HELLS", "CCNB1",
                   "E2F2", "CDC45", "OIP5", "GINS1", "MCM3", "MCM4", "MCM5", "RECQL4", "PCNA", "KIF20A", "RFC3", 
                   "KPNA2", "NCAPH", "TONSL", "ZWILCH", "CHAF1A", "MKI67", "ESPL1", "FANCD2", "FEN1", "FOXM1",
                   "CCl1", "LTA", "XCL1", "XCL2", "CCL4", "CSF2", "IFNG", "TNF", "CCL4L2", "PRF1", "GZMB", "GZMH",
                   "VSIR", "CD96", "TIGIT", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "CRTAM", 
                   "CD226", "CD244", "KIR2DL4", "NCR1", "NCR2", "NCR3")
  x <- LGGno0 %>%
    spread(celltype, scale) %>%
    dplyr::select(sample, contains("nk_"), contains("t_"), -mast_cell) %>%
    right_join(foreach(i = length(genelist), .combine = bind_rows) %do%{
      if(i < 42){n <- LGG %>% 
        filter(symbol == genelist[[i]]) %>% mutate(cat = "Cell Cycle", item = log(raw_count_scaled + 1))}
      else if(i < 51 & i >= 42){n <- LGG %>% 
        filter(symbol == genelist[[i]]) %>% mutate(cat = "Chemokine", item = log(raw_count_scaled + 1))}
      else if(i >= 51 & i < 54){n <- LGG %>% 
        filter(symbol == genelist[[i]]) %>% mutate(cat = "Cytotoxicity", item = log(raw_count_scaled + 1))}
      else if(i >= 54 & i < 61){n <- LGG %>% 
        filter(symbol == genelist[[i]]) %>% mutate(cat = "Check Point", item = log(raw_count_scaled + 1))}
      else if(i >= 61){n <- LGG %>% 
        filter(symbol == genelist[[i]]) %>% mutate(cat = "Receptor", item = log(raw_count_scaled + 1))}
      n
    } %>%
      dplyr::select(sample, symbol, cat, item) %>%
      spread(symbol, item))
  
  
  a <- melt(cor(x %>% dplyr::select(-sample))[-c(1:8), c(1:8)])
  a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
  a <- a %>% inner_join(as.data.frame(a.p$P[-c(1:8), c(1:8)]) %>% 
                          mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:8), c(1:8)]))) %>% 
                          gather(X2, p.value, -X1)) 
  
  a$stars <- cut(a$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
  a$X2 <- factor(a$X2, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD', "t_helper", "t_CD8_naive", "t_gamma_delta", 't_CD4_memory', "t_CD8_memory"), ordered=TRUE)
  
  ggplot(data = a , aes(x=`X2`, y=`X1`, fill=value)) + 
    geom_tile() + 
    geom_text(aes(label=stars)) + 
    scale_fill_viridis() +
    labs(y = "", x = "") + #y = "Log Transformed Normalized Count", x = "Log Transformed NK Proportion"
    theme_classic() +
    scale_x_discrete(labels = c("ReNK", "IL2NK", "SPANK", "Helper T", "Naive CD8 T", "?æ? T", "Memory CD4 T", "Memory CD8 T")) +
    theme(axis.text.x = element_text(angle = -45)) +
    labs(fill = "Correlation")
}
