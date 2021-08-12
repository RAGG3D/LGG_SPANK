#Fig 7
CGGA <- data.table::fread("data/CGGA/mRNA_1/CGGA.mRNAseq_325.RSEM-genes.20200506.txt") %>%
  pivot_longer(contains("CGGA"), names_to = "sample", values_to = "RSEM") %>%
  tidybulk::rename("symbol" = "Gene_Name") %>%
  inner_join(data.table::fread("data/CGGA/mRNA_1/CGGA.mRNAseq_325_clinical.20200506.txt") %>%
               tidybulk::rename("sample" = "CGGA_ID")) %>%
  bind_rows(data.table::fread("data/CGGA/mRNA_2/CGGA.mRNAseq_693.RSEM-genes.20200506.txt") %>%
              pivot_longer(contains("CGGA"), names_to = "sample", values_to = "RSEM") %>%
              tidybulk::rename("symbol" = "Gene_Name") %>%
              inner_join(data.table::fread("data/CGGA/mRNA_2/CGGA.mRNAseq_693_clinical.20200506.txt") %>%
                           tidybulk::rename("sample" = "CGGA_ID"))) %>%
  filter(Grade == "WHO II") %>%
  tidybulk::rename("vital_status" = "Censor (alive=0; dead=1)", "total_living_days" = "OS") %>%
  na.omit()

CGGA_cell <- read_csv("data/CGGA/CGGA_old_cibersort.csv") %>%
  tidybulk::rename(sample = `Input Sample`) %>%
  dplyr::select(-`P-value`, -`Pearson Correlation`, -RMSE) %>%
  tidybulk::rename("ReNK" = "nk_resting",
                   "IL2NK" = "nk_primed_IL2",
                   "SPANK" = "nk_primed_IL2_PDGFD",
                   "Memory CD8 T" = "t_CD8_memory",
                   "Helper T" = "t_helper")

######FIG A
x <- foreach(i = list("PDGFD", "TGFBI", "IGFBP3", "CHI3L1"), .combine = bind_rows) %do%{
  CGGA %>% 
    filter(symbol == i) %>%
    mutate(item = log(RSEM + 0.1), cat = symbol) 
} %>%
  dplyr::select(sample, cat, item) %>%
  spread(cat, item)
a <- melt(cor(x %>% dplyr::select(-sample))[-c(1,2,4), c(1,2,4)])
a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
a <- a %>% 
  mutate(X1 = rownames(.)) %>%
  inner_join(as.data.frame(a.p$P[-c(1,2,4), c(1,2,4)]) %>% 
               mutate(X1 = rownames(as.data.frame(a.p$P[-c(1,2,4), c(1,2,4)]))) %>% 
               gather(X2, p.value, -X1)) 

a$stars <- cut(a$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
a <- a %>% mutate(X2 = "PDGFD")

pa01 <- ggplot(a, aes(x = X1, y = log(p.value*10000 + 0.00001), group_by(X1))) +
  geom_point(size = 3) +
  facet_wrap(.~X2, nrow = 2) +
  labs(x = "", y = "Log(pval*10000)") +
  geom_hline(yintercept=log(0.05*10000 + 0.00001), linetype="dashed") +
  theme_bw() + 
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) 


x <- as.data.frame( CGGA %>% 
                      filter(symbol == "PDGFD") %>%
                      dplyr::select(sample, symbol, RSEM) %>%
                      spread(symbol, RSEM) %>%
                      mutate(PDGFD = factor(Hmisc::cut2(PDGFD, g = 2), labels = c("L", "H"))) %>%
                      inner_join(data.table::fread("data/CGGA/mRNA_2/CGGA.mRNAseq_693_clinical.20200506.txt") %>%
                                   bind_rows(data.table::fread("data/CGGA/mRNA_1/CGGA.mRNAseq_325_clinical.20200506.txt")) %>%
                                   tidybulk::rename("sample" = "CGGA_ID"))%>%
                      filter(Grade == "WHO II") %>%
                      tidybulk::rename("vital_status" = "Censor (alive=0; dead=1)", "total_living_days" = "OS") %>%
                      na.omit()) %>%
  mutate(cat = "PDGFD")

x$PDGFD <- factor(x$PDGFD, levels=c("L", "H"), ordered=TRUE)

pa02 <- ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ PDGFD,
    data = x 
  ),
  data = x,
  facet.by = c("cat"),
  conf.int = T,
  risk.table = F,
  conf.int.alpha = 0.15,
  pval = T,
  linetype = "PDGFD",
  legend.title = "Abundance ",
  legend.labs = list("L ", "H "),
  short.panel.labs = T,
  palette = "jco"
) + theme_bw() + 
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE)


######FIG B
x <-as.data.frame(   
  foreach(i = c("ReNK", "IL2NK", "SPANK"), .combine = bind_rows) %do%{
    
    CGGA %>%
      filter(symbol == "PDGFD") %>%
      dplyr::select(sample, symbol, RSEM) %>%
      spread(symbol, RSEM) %>%
      mutate(gene = factor(Hmisc::cut2(PDGFD, g = 2), labels = c(1:nlevels(Hmisc::cut2(PDGFD, g = 2)))))  %>%
      inner_join(CGGA_cell %>%
                   mutate(cell = factor(Hmisc::cut2(!!as.name(i), g = 2), labels = c(1:nlevels(Hmisc::cut2(!!as.name(i), g = 2))))) %>%
                   dplyr::select(sample, cell)) %>%
      unite("item", gene, cell, sep = "/", remove = T) %>%
      dplyr::select(sample, item) %>% 
      mutate(cat = paste0("PDGFD/", i)) %>%
      mutate(item = gsub("1", "L", item),
             item = gsub("2", "H", item))
  }) %>%
  inner_join(data.table::fread("data/CGGA/mRNA_2/CGGA.mRNAseq_693_clinical.20200506.txt") %>%
               bind_rows(data.table::fread("data/CGGA/mRNA_1/CGGA.mRNAseq_325_clinical.20200506.txt")) %>%
               tidybulk::rename("sample" = "CGGA_ID"))%>%
  filter(Grade == "WHO II") %>%
  tidybulk::rename("vital_status" = "Censor (alive=0; dead=1)", "total_living_days" = "OS") %>%
  na.omit()

x$item <- factor(x$item, levels=c("L/L", "L/H", "H/L", "H/H"), ordered=TRUE)
x$cat <- factor(x$cat, levels=c("PDGFD/ReNK", "PDGFD/IL2NK", "PDGFD/SPANK"), ordered=TRUE)

pb <- ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("cat"),
  conf.int = T,
  legend.title = "Gene Abundance/NK",
  legend.labs = list("L/L", "L/H", "H/L", "H/H"),
  conf.int.alpha = 0.15,
  risk.table = F,
  pval = T,
  short.panel.labs = T,
  linetype = "item",
  palette = "jco"
) + theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) 


######FIG C
x <- as.data.frame( foreach(i = list("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4"), .combine = bind_rows) %do% {
  CGGA %>% filter(symbol == i) %>%
    mutate(median = factor(Hmisc::cut2(RSEM, g = 2), labels = c(1:nlevels(Hmisc::cut2(RSEM, g = 2))))) 
} %>% 
  dplyr::select(sample, symbol, median) %>%
  gather(cat, item, -c(sample, symbol)) %>%
  inner_join(data.table::fread("data/CGGA/mRNA_2/CGGA.mRNAseq_693_clinical.20200506.txt") %>%
               bind_rows(data.table::fread("data/CGGA/mRNA_1/CGGA.mRNAseq_325_clinical.20200506.txt")) %>%
               tidybulk::rename("sample" = "CGGA_ID"))%>%
  filter(Grade == "WHO II") %>%
  tidybulk::rename("vital_status" = "Censor (alive=0; dead=1)", "total_living_days" = "OS") %>%
  na.omit() %>%
  mutate(item = gsub("1", "L", item),
         item = gsub("2", "H", item)))

x$symbol <- factor(x$symbol, levels=c("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4"), ordered=TRUE)
x$item <- factor(x$item, levels=c("L", "H"), ordered=TRUE)
pc <- ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("symbol"),
  conf.int = T,
  risk.table = F,
  short.panel.labs = T,
  pval = T,
  nrow =1,
  conf.int.alpha = 0.15,
  legend.title = "Abundance",
  palette = "jco"
) + theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") 

######FIG D
cell <- CGGA_cell %>%
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

Immunereceptor <- list("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4")
x <- no0cell %>%
  spread(celltype, scale) %>%
  dplyr::select(sample, contains("NK"), `Helper T`, `Memory CD8 T`, -mast_cell) %>%
  inner_join(foreach(i = Immunereceptor, .combine = bind_rows) %do%{
    CGGA %>% 
      filter(symbol == i) %>%
      mutate(item = log(RSEM + 1), cat = symbol) 
  } %>%
    dplyr::select(sample, cat, item) %>%
    spread(cat, item))

a <- melt(cor(x %>% dplyr::select(-sample))[-c(1:5), c(1:5)])
a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
pd <- as_tibble(a %>% inner_join(as.data.frame(a.p$P[-c(1:5), c(1:5)]) %>% 
                                   mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:5), c(1:5)]))) %>% 
                                   gather(X2, p.value, -X1)) ) %>%
  tidybulk::rename("Cell Type" = "X2",
                   "Receptor" = "X1") %>%
  tidyHeatmap::heatmap(
    Receptor,
    `Cell Type`,
    value,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    palette_value = circlize::colorRamp2(
      
      seq(0.5,-0.5, length.out = 11),
      
      RColorBrewer::brewer.pal(11, "RdBu")
      
    )
    
  ) 


ggsave(plot = plot_grid(plot_grid(pa01, pa02, pb, nrow = 1, align = "hv", axis = "b", rel_widths = c(1.4,1.4,3.4)), 
                        pc, ncol = 1, align = "hv", axis = "t"), "output/LGG-p6.pdf", height = 10, width = 18)
