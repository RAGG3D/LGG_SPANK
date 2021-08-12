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
