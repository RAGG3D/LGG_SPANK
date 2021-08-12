#FIG 5
Tcell <- c("t_CD4_memory",  "t_CD8_memory",  "t_CD8_naive",   "t_gamma_delta", "t_helper")
x <- foreach(i = 1:length(Tcell), .combine = bind_rows) %do% {
  read_csv(paste0("data/TCGA_LGG_CIBERSORT_PDGFDD.csv")) %>%
    tidybulk::rename(sample = `Input Sample`) %>%
    dplyr::select(-`P-value`, -`Pearson Correlation`, -RMSE) %>%
    gather(cat, item, -sample) %>% 
    filter(cat == Tcell[i]) %>%
    mutate(item = ifelse(item > median(item), "H", "L"))} %>%
  TCGA_clinical("LGG")
x <- as.data.frame(x)
x$item <- factor(x$item, levels=c("L", "H"), ordered=TRUE)
x$cat <- factor(x$cat, levels=c("t_helper", "t_CD8_naive", "t_gamma_delta", 't_CD4_memory', "t_CD8_memory"), ordered=TRUE)
p5a =  ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("cat"),
  conf.int = T,
  conf.int.alpha = 0.15,
  risk.table = F,
  nrow = 1,
  pval = F,
  legend.title = "T cell Infiltration Fraction",
  legend.labs = list("L ", "H "),
  short.panel.labs = T,
  panel.labs = list(cat = c("Helper T", "Naive CD8 T", "?æ? T", "Memory CD4 T", "Memory CD8 T")),
  palette = "jco"
) + theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom")

x <- Gene_Cell("LGG", "PDGFD", c("t_helper", "t_CD8_naive", "t_gamma_delta", 't_CD4_memory', "t_CD8_memory")) %>%
  rbind(Gene_Cell("LGG", "PDGFRB", c("t_helper", "t_CD8_naive", "t_gamma_delta", 't_CD4_memory', "t_CD8_memory"))) %>%
  mutate(cat = gsub("t_helper", "Helper T", cat)) %>%
  mutate(cat = gsub("t_CD8_naive", "Naive CD8 T", cat)) %>%
  mutate(cat = gsub("t_gamma_delta", "GD T", cat)) %>%
  mutate(cat = gsub("t_CD4_memory", "Memory CD4 T", cat)) %>%
  mutate(cat = gsub("t_CD8_memory", "Memory CD8 T", cat)) 
x <- as.data.frame(x)
x <- x %>% 
  #filter(total_living_days >= 0) %>%
  filter( grepl('PDGFD/', cat))
x$item <- factor(x$item, levels=c("1/1", "1/2", "2/1", "2/2"), ordered=TRUE)
x$cat <- factor(x$cat, levels=c("PDGFD/Helper T", "PDGFD/Naive CD8 T", "PDGFD/GD T", "PDGFD/Memory CD4 T", "PDGFD/Memory CD8 T"), ordered=TRUE)

p5b =  ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("cat"),
  conf.int = T,
  nrow = 1,
  legend.title = "Gene Abundance/T",
  legend.labs = list("L/L", "L/H", "H/L", "H/H"),
  conf.int.alpha = 0.15,
  risk.table = F,
  pval = F,
  short.panel.labs = T,
  linetype = "item",
  palette = "jco"
) + theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) 

p6 <- SPANK_T_cor("LGG") + theme(text=element_text(size=16, family="sans"),
                                 legend.position = "bottom") 

ggsave(plot = plot_grid(p5a, p5b,p6,
                        ncol = 1), "output/LGG-p5.pdf",device = "pdf", height = 12, width = 15)

