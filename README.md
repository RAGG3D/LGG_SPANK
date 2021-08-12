# LGG_SPANK

## Packages, functions and data required

```r
source("src/functions.R") 

```

## Figure 1 (analysis from @stemangiola)

All the figure should be saved as a PDF for further edits.
(After modification)

![image](https://github.com/RAGG3D/LGG_SPANK/blob/main/figure/LGG-p1%2020cm.jpg)

## Figure 2 (See Figure scripts/Figure 2.R)
### FIG 2A
```r
x <- Bar_TCGAcelltype("LGG")
x$nkstate = factor(x$nkstate, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD'))
p2a01 <- x %>%
  filter(cat == "Fraction") %>%
  ggplot(aes(x = reorder(order, x = sample), y = profile, fill = nkstate)) +
  geom_bar(stat = "identity") +
  #facet_grid(~cat, scales = "free") +
  scale_fill_brewer(palette="Blues") +
  coord_flip() + 
  theme_bw() +
  theme(panel.background=element_rect(fill='transparent',color ="gray")) +
  labs( x = "Patients", y = "Fraction") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text=element_text(size=16, family="sans"),
        legend.position = "bottom")
p2a02 <- x %>%
  filter(cat == "Percentage") %>%
  ggplot(aes(x = reorder(order, x = sample), y = profile, fill = nkstate)) +
  geom_bar(stat = "identity") +
  #facet_grid(~cat, scales = "free") +
  scale_fill_brewer(palette="Blues", labels = c("ReNK", "IL2NK", "SPANK")) +
  coord_flip() + 
  theme_bw() +
  theme(panel.background=element_rect(fill='transparent',color ="gray")) +
  labs(x = "Patients", y = "Percentage") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        text=element_text(size=16, family="sans"),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") + 
  guides(fill=guide_legend(title="NK Phenotype"))
 ```

### FIG 2B
```r
x <- KM_TCGAcelltype_median("LGG")
x <- as.data.frame(x)
x$celltype = factor(x$celltype, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD'))
x$fraction = factor(x$fraction, levels=c("L", "H"))
p2b =  ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ fraction,
    data = x 
  ),
  data = x,
  facet.by = c("celltype"),
  conf.int = T,
  risk.table = F,
  conf.int.alpha = 0.15,
  pval = F,
  legend.title = "Cell Fraction",
  legend.labs = list("L ", "H "),
  short.panel.labs = T,
  panel.labs = list(celltype = c("ReNK", "IL2NK", "SPANK")),
  linetype = "fraction",
  palette = "jco"
) + theme_bw() + 
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) 
```
Then we need to format the 2 panels:
```r
plot_grid(plot_grid(p2a01, p2a02, nrow = 1), p2b, ncol = 1)
```
(After modification)

![image](https://github.com/RAGG3D/LGG_SPANK/blob/main/figure/LGG-Fig2.jpg)

## Figure 3 (See Figure scripts/Figure 3.R)
### FIG 3A
```r
gene <- list("PDGFD", "PDGFRB", "TGFBI", "IGFBP3", "CHI3L1")
p3a <- Gene_heatmap_log10("LGG", gene) 
x <- Gene_correlation("LGG", gene)
p3a <- ggplot(x, aes(x = X1, y = log(p.value*10000 + 0.00001), group_by(X1))) +
  geom_point(size = 3) +
  facet_wrap(.~X2, nrow = 2) +
  labs(x = "", y = "Log(pval*10000)") +
  geom_hline(yintercept=log(0.05*10000 + 0.00001), linetype="dashed") +
  theme_bw() + 
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) 
```

### FIG 3B
```r
x <- PDGFsurvival("LGG", 2) %>%
  gather(cat, item, -sample) %>% 
  TCGA_clinical("LGG")
x <- as.data.frame(x)
x$cat <- factor(x$cat, levels=c("PDGFD", "PDGFRB"), ordered=TRUE)
x$item <- factor(x$item, levels=c(1, 2), ordered=TRUE)
p3b =  ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("cat"),
  conf.int = T,
  risk.table = F,
  conf.int.alpha = 0.15,
  pval = F,
  linetype = "item",
  legend.title = "Abundance ",
  legend.labs = list("L ", "H "),
  short.panel.labs = T,
  palette = "jco"
) + theme_bw() + 
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) +
  theme(aspect.ratio=1)
```

### FIG 3C
```r
x <- PDGFsurvival("LGG", 2) %>%
  unite("item", PDGFD:PDGFRB, remove = T, sep = "/") %>%
  mutate(cat = "PDGFD/PDGFRB") %>%
  TCGA_clinical("LGG")
x <- as.data.frame(x)
x$item <- factor(x$item, levels=c("1/1", "1/2", "2/1", "2/2"), ordered=TRUE)
p3c =  ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("cat"),
  conf.int = T,
  risk.table = F,
  conf.int.alpha = 0.15,
  short.panel.labs = T,
  pval = F,
  linetype = "item",
  legend.title = "",
  legend.labs = list("L/L", "L/H", "H/L", "H/H"),
  palette = "jco"
) + theme_bw() + 
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) +
  theme(aspect.ratio=1)
```
Then we need to format the 3 panels:
```r
plot_grid(p3a, p3b, p3c, nrow = 1, rel_widths = c(1.5,2.5,1.5))
```

(After modification)

![image](https://github.com/RAGG3D/LGG_SPANK/blob/main/figure/LGG-p3%2020cm.jpg)

## Figure 4 (See Figure scripts/Figure 4.R)
### FIG  4
```r
x <- Gene_Cell("LGG", "PDGFD", c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD')) %>%
  rbind(Gene_Cell("LGG", "PDGFRB", c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD'))) %>%
  mutate(cat = gsub("nk_resting", "ReNK", cat)) %>%
  mutate(cat = gsub("nk_primed_IL2", "IL2NK", cat)) %>%
  mutate(cat = gsub("IL2NK_PDGFD", "SPANK", cat)) 
x <- as.data.frame(x)
x$item <- factor(x$item, levels=c("1/1", "1/2", "2/1", "2/2"), ordered=TRUE)
x <- x %>% 
  #filter(total_living_days >= 0) %>%
  filter( grepl('PDGFD/', cat))
x$cat <- factor(x$cat, levels=c("PDGFD/ReNK", "PDGFD/IL2NK", "PDGFD/SPANK"), ordered=TRUE)



p4 =  ggsurvplot(
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
  pval = F,
  short.panel.labs = T,
  linetype = "item",
  palette = "jco"
) + theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") +
  guides(linetype = FALSE) 
```
(After modification)

![image](https://github.com/RAGG3D/LGG_SPANK/blob/main/figure/LGG-p4%2020cm.jpg)

## Figure 5 (See Figure scripts/Figure 5.R)
### FIG 5A
```r
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
```

### FIG 5B
```r
x <- Gene_Cell("LGG", "PDGFD", c("t_helper", "t_CD8_naive", "t_gamma_delta", 't_CD4_memory', "t_CD8_memory")) %>%
  rbind(Gene_Cell("LGG", "PDGFRB", c("t_helper", "t_CD8_naive", "t_gamma_delta", 't_CD4_memory', "t_CD8_memory"))) %>%
  mutate(cat = gsub("t_helper", "Helper T", cat)) %>%
  mutate(cat = gsub("t_CD8_naive", "Naive CD8 T", cat)) %>%
  mutate(cat = gsub("t_gamma_delta", "GD T", cat)) %>%
  mutate(cat = gsub("t_CD4_memory", "Memory CD4 T", cat)) %>%
  mutate(cat = gsub("t_CD8_memory", "Memory CD8 T", cat)) 
x <- as.data.frame(x)
x <- x %>% 
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
```

### FIG 5C
```r
p5c <- SPANK_T_cor("LGG") + theme(text=element_text(size=16, family="sans"),
                                 legend.position = "bottom") 
```

Then we need to format the panels:
```r
plot_grid(p5a, p5b, p5c, ncol = 1)
```
and modify the figure:

(After modification)

![image](https://github.com/RAGG3D/LGG_SPANK/blob/main/figure/LGG-p5%2020cm.jpg)

## Figure 6 (See Figure scripts/Figure 6.R)
### FIG 6A
```r
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
```

### FIG 6B
```r
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

```
We need to save p6b seperately, because it a InputHeatmap instead of a grob. So p5c cannot be formatted by plot_grid with p6a, and modify the figure:

(After modification)

![image](https://github.com/RAGG3D/LGG_SPANK/blob/main/figure/LGG-p6%2020cm.jpg)

## Figure 7 (See Figure scripts/Figure 7.R)
### Figure 7 is based on another dataset, CGGA, so we need to format the data matrix firstly:  
```r
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
```

### FIG 7A
```r
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
```

### FIG 7B
```r
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
```

### FIG 7C
```r
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
```

### FIG 7D
```r
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
```
We need to save p7d seperately, because it a InputHeatmap instead of a grob. So p5c cannot be formatted by plot_grid with other panels, and modify the figure:

(After modification)

![image](https://github.com/RAGG3D/LGG_SPANK/blob/main/figure/LGG-p7%2020cm.jpg)


