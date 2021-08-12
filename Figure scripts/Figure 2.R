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
