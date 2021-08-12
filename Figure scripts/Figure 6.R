#FIG 6
x <- genelist_cancer(LGGBLCA[["lgg"]], "LGG", list("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "NCR2")) %>%
  filter(cat == "median")
x$symbol <- factor(x$symbol, levels=c("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "NCR2"), ordered=TRUE)
p6a <- ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("symbol"),
  conf.int = T,
  risk.table = F,
  legend.labs = list("L", "H", "?", "?"),
  short.panel.labs = T,
  pval = F,
  nrow =2,
  conf.int.alpha = 0.15,
  legend.title = "Abundance",
  palette = "jco"
) + theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") 

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

a <- melt(cor(x %>% dplyr::select(-sample))[-c(1:5), c(1:5)]) %>%
  mutate(X2 = gsub("nk_resting", "ReNK", X2)) %>%
  mutate(X2 = gsub("nk_primed_IL2_PDGFD", "SPANK", X2)) %>%
  mutate(X2 = gsub("nk_primed_IL2", "IL2NK", X2)) %>%
  mutate(X2 = gsub("t_CD8_memory", "Memory CD8 T", X2)) %>%
  mutate(X2 = gsub("t_helper", "Helper T", X2)) %>%
  mutate(X2 = fct_relevel(X2, "ReNK", "IL2NK", "SPANK", "Helper T", "Memory CD8 T")) %>%
  mutate(X1 = fct_relevel(X1, "NCR3", "NCR2","NCR1","KLRK1",  "KLRC4","KLRC3","KLRC2","KLRC1", "KIR2DL4", "CRTAM",  "CD244",
                          "CD226"))

p6b <- a %>%
  as_tibble() %>%
  tidybulk::rename("Cell Type" = "X2",
                   "Receptor" = "X1") %>%
  tidyHeatmap::heatmap(
    Receptor,
    `Cell Type`,
    value,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    palette_value = circlize::colorRamp2(
      
      seq(2,-2, length.out = 11),
      
      RColorBrewer::brewer.pal(11, "RdBu")
      
    )
    
  ) %>% tidyHeatmap::save_pdf("output/LGG-p6-heatmap.pdf", width=4, height = 10)
