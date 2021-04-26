options(connectionObserver = NULL)

library(tidyverse)
library(tidybulk)
library(survminer)
library(survival)
library(gapminder)
library(foreach)
library(org.Hs.eg.db)
library(cowplot)
library(ggsci)
library(GGally)
library(gridExtra)
library(grid)
library(reshape)
library(Hmisc)
library(viridis)
library(furrr)

source("src/functions.R")
LGG <- LGG_transcript() %>%
filter(symbol == "PDGFD"|symbol == "PDGFRB"|symbol == "KLRK1"|symbol == "KLRC1"|
         symbol == "KLRC2"|symbol == "KLRC3"|symbol == "KLRC4"|symbol == "NCR2"|
         symbol == "NCR1"|symbol == "NCR3"|symbol == "KIR2DL4"|symbol == "CRTAM"|
         symbol == "CD244"|symbol == "CD226"|symbol == "TGFBI"|symbol == "IGFBP3"|
         symbol == "CHI3L1") %>%
  clinical_combine("LGG") %>% 
  dplyr::select(sample, symbol, raw_count_scaled, vital_status, total_living_days, age)

LGGno0 <- no0_CIBERSORT("LGG")

#FIG 2
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

ggsave(plot = plot_grid(plot_grid(p2a01, p2a02, nrow = 1), p2b, ncol = 1), "output/LGG-p2.pdf", device = "pdf", height = 7, width = 8)


#FIG 3
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

ggsave(plot = plot_grid(p3a, p3b, p3c, nrow = 1, rel_widths = c(1.5,2.5,1.5)),
       "output/LGG-p3.pdf", device = "pdf", height = 4, width = 14)

x <- foreach(i = list_TCGA, .combine = bind_rows) %do% {
  x <- no0_CIBERSORT(i) %>%
    spread(celltype, scale) %>%
    dplyr::select(sample, contains("nk")) 
  
  a <- melt(cor((x)[, -c(1)]))
  a.p <- rcorr(as.matrix(x)[, -c(1)])
  a <- a %>% inner_join(as.data.frame(a.p$P) %>% 
                          mutate(X1 = rownames(as.data.frame(a.p$P))) %>% 
                          gather(X2, p.value, -X1)) %>%
    tidybulk::rename(Correlation = value) %>%
    mutate(cancer = i)
  a
}
a <- x %>%
  mutate(X1 = gsub("nk_resting", "ReNK", X1)) %>%
  mutate(X1 = gsub("nk_primed_IL2_PDGFD", "SPANK", X1)) %>%
  mutate(X1 = gsub("nk_primed_IL2", "IL2NK", X1)) %>%
  mutate(X2 = gsub("nk_resting", "ReNK", X2)) %>%
  mutate(X2 = gsub("nk_primed_IL2_PDGFD", "SPANK", X2)) %>%
  mutate(X2 = gsub("nk_primed_IL2", "IL2NK", X2)) %>%
  arrange(cancer)

a$stars <- cut(a$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
a$X2 = factor(a$X2, levels=c('ReNK', 'IL2NK', 'SPANK'))
a$X1 = factor(a$X1, levels=c('ReNK', 'IL2NK', 'SPANK'))

p3d01 <- ggplot(data = a[c(1:153),] , aes(x=X2, y=X1, fill=Correlation)) + 
  geom_tile() + 
  geom_text(aes(label=stars)) + 
  scale_fill_viridis() +
  #scale_fill_gradient2(high="#FC4E07", mid="#E7B800", low="#00AFBB") +  
  #scale_fill_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) +
  labs(y = "", x = "") + #y = "Log Transformed Normalized Count", x = "Log Transformed NK Proportion"
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(.~cancer, space = "free_y", scales="free_y", switch="y") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill=NA,colour=NA),
        panel.spacing=unit(0,"cm"), axis.title.y = element_blank(),
        strip.text.y = element_text(size = 11)) +
  coord_cartesian(clip = "off")

ggsave(plot = p3d01, "G:/Thesis/LGG paper plot/LGG-p3d01.png", height = 3, width = 15)

#FIG 4
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

ggsave(plot = p4, "output/LGG-p4.pdf",device = "pdf", height = 4, width = 10)

#PDGFD/SPANK
x <- x %>% 
  filter(cat == "PDGFD/SPANK") %>%
  filter(item == "1/1"|item == "1/2") %>%
  mutate(subcat = "PDGFD Low") %>%
  rbind(x %>% 
          filter(cat == "PDGFD/SPANK") %>%
          filter(item == "2/1"|item == "2/2") %>%
          mutate(subcat = "PDGFD High")) %>%
  rbind(x %>% 
          filter(cat == "PDGFD/SPANK") %>%
          mutate(item = ifelse(item == "1/2", item, "others")) %>%
          mutate(subcat = "L/H to others")) %>%
  rbind(x %>% 
          filter(cat == "PDGFD/SPANK") %>%
          mutate(subcat = "Combination"))

x <- as.data.frame(x)
x$subcat <- factor(x$subcat, levels=c("Combination", "PDGFD Low", "PDGFD High", "L/H to others"), ordered=TRUE)
x$item <- factor(x$item, levels=c("1/1", "1/2", "2/1", "2/2", "others"), ordered=TRUE)


p =  ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("cat", "subcat"),
  conf.int = T,
  legend.title = "Gene Abundance/NK",
  legend.labs = list("L/L", "L/H", "H/L", "H/H", "others"),
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



ggsave(plot = p4, "G:/Thesis/LGG paper plot/LGG-p4.png", height = 7, width = 10)

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


#FIG 7
x <- genelist_cancer(LGGBLCA[["lgg"]], "LGG", list("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "NCR2")) %>%
  filter(cat == "median")
x$symbol <- factor(x$symbol, levels=c("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "NCR2"), ordered=TRUE)
p7a <- ggsurvplot(
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

x <- genelist_cancer(LGGBLCA[["lgg"]], "LGG", list("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "NCR2")) %>%
  filter(cat == "quantile")
x$symbol <- factor(x$symbol, levels=c("KLRK1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "NCR2"), ordered=TRUE)
p7b <- ggsurvplot(
  survfit(
    Surv(total_living_days, vital_status) ~ item,
    data = x 
  ),
  data = x,
  facet.by = c("symbol"),
  conf.int = T,
  risk.table = F,
  legend.labs = list("Q1", "Q2", "Q3", "Q4"),
  short.panel.labs = T,
  pval = F,
  nrow =2,
  conf.int.alpha = 0.15,
  legend.title = "Abundance",
  palette = "jco"
) + theme_bw() +
  labs(tag = "B") +
  theme(text=element_text(size=16, family="sans"),
        legend.position = "bottom") 

p7c <- NK_selectT_receptor_heatmap_log10("LGG") +
  labs(tag = "C")
ggsave(plot = plot_grid(p7a, p7b, p7c, ncol = 1), "G:/Thesis/LGG paper plot/LGG-p7.pdf",device = "pdf", height = 19, width = 9)


#FIG Supplementary 1
y <- foreach(i = c("t_CD4_memory",  "t_CD8_memory",  "t_CD8_naive",   "t_gamma_delta", "t_helper",
                   "nk_resting", "nk_primed_IL2", "nk_primed_IL2_PDGFD"), .combine = inner_join) %do%{
                     read_csv("data/TCGA_LGG_CIBERSORT_PDGFDD.csv") %>%
                       mutate(sample = `Input Sample`) %>%
                       dplyr::select(sample, contains(i))
                   } %>%
  gather(celltype, fraction, -sample) 
x <- y %>%
  group_by(sample) %>% 
  summarise(order = sum(fraction)) %>%
  right_join(y, by = c("sample"))

x <- as.data.frame(x)
x$celltype <- factor(x$celltype, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD', "t_helper", "t_CD8_naive", "t_gamma_delta", 't_CD4_memory', "t_CD8_memory"), ordered=TRUE)

ps1a <- x %>%
  ggplot(aes(fill=celltype, y=fraction, x=reorder(order, x = sample))) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(name = "Cell Type", labels = c("ReNK", "IL2NK", "SPANK", "Helper T", "Naive CD8 T", "?æ? T", "Memory CD4 T", "Memory CD8 T"), 
                    values = c("#2c9321", "#a4db77", "#1b63a5", "#96c3dc", "#d90017", "#f88587", "#fc6908", "#fbb25c")) +
  theme_classic() +
  labs(tag = "A", y = "Fraction", x = "Patients") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 

ps1b <- x %>%
  ggplot(aes(fill=celltype, y=fraction, x=reorder(order, x = sample))) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(name = "Cell Type", labels = c("ReNK", "IL2NK", "SPANK", "Helper T", "Naive CD8 T", "?æ? T", "Memory CD4 T", "Memory CD8 T"), 
                    values = c("#2c9321", "#a4db77", "#1b63a5", "#96c3dc", "#d90017", "#f88587", "#fc6908", "#fbb25c")) +
  theme_classic() +
  #theme(panel.background=element_rect(fill='transparent',color ="gray")) +
  labs(tag = "B", x = "Patients", y = "Percentage") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 

ggsave(plot = plot_grid(ps1a, ps1b, ncol = 1), "output/LGG-FigSup1.png", height = 6, width = 8)

#FIG Supplementary 2
x <- genelist_cancer(LGG, "LGG", list("CD226", "CD244", "CRTAM", "KIR2DL4", "NCR1", "NCR3")) %>%
  filter(cat == "median")
x$symbol <- factor(x$symbol, levels=c("CD226", "CD244", "CRTAM", "KIR2DL4", "NCR1", "NCR3"), ordered=TRUE)
ps2 <- ggsurvplot(
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
) + theme_bw() 

ggsave(plot = ps2, "output/LGG-FigSup2.png", height = 6, width = 10)

#Supplementary material 4
foreach(i = list_TCGA) %do% {
  x <- Bar_TCGAcelltype(i)
  x$nkstate = factor(x$nkstate, levels=c('nk_resting', 'nk_primed_IL2', 'nk_primed_IL2_PDGFD'))
  p2a01 <- x %>%
    filter(cat == "Fraction") %>%
    ggplot(aes(x = reorder(order, x = sample), y = profile, fill = nkstate)) +
    geom_bar(stat = "identity") +
    #facet_grid(~cat, scales = "free") +
    scale_fill_brewer(palette="Blues") +
    coord_flip() + 
    theme_classic() +
    theme(panel.background=element_rect(fill='transparent',color ="gray")) +
    labs(tag = "A", x = "Patients", y = "Fraction", title = paste0("TCGA-", i)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
  p2a02 <- x %>%
    filter(cat == "Percentage") %>%
    ggplot(aes(x = reorder(order, x = sample), y = profile, fill = nkstate)) +
    geom_bar(stat = "identity") +
    #facet_grid(~cat, scales = "free") +
    scale_fill_brewer(palette="Blues", labels = c("ReNK", "IL2NK", "SPANK")) +
    coord_flip() + 
    theme_classic() +
    theme(panel.background=element_rect(fill='transparent',color ="gray")) +
    labs(tag = "", x = "Patients", y = "Percentage", title = "") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) + 
    guides(fill=guide_legend(title="NK Phenotype"))
  
  x <- KM_TCGAcelltype_median(i)
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
    guides(linetype = FALSE) + 
    labs(tag = "B", title = "")
  
  ggsave(plot = plot_grid(plot_grid(p2a01, p2a02, rel_widths = c(0.67, 1)), p2b, nrow = 1,
                          rel_widths = c(1.2,2)), paste0("output/", i, ".png"), height = 3.3, width = 14)
  
}