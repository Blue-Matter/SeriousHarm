
SH <- readxl::read_excel('Report/Seriousharm_estimates_import.xlsx') %>% select(!`F/F[rep]`) %>%
  filter(Stock != "US Pacific ocean perch" & Stock != "GM haddock")

SH <- SH[order(SH$`B[ESH]/B[init]`), ] %>% 
  #mutate(St = paste0("(", 1:nrow(SH), ")")) %>% 
  mutate(St = 1:nrow(SH)) 

SH_order <- SH %>% 
  reshape2::melt(id.vars = c("Stock", "St", "Year_assess", "B[ESH]/B[init]", "tv")) %>%
  mutate(value = ifelse(value > 5, 5, value)) %>%
  filter(!is.na(value))

g <- ggplot(SH_order, aes(`B[ESH]/B[init]`, value)) + 
  #geom_smooth() + 
  geom_text(aes(label = St), size = 3) + 
  #geom_point() +
  #ggrepel::geom_text_repel(aes(label = St)) + 
  theme_bw() + 
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA, colour = NA), 
        axis.title.y = element_blank(),
        strip.placement = "outside") + 
  facet_wrap(vars(variable), scales = "free_y", labeller = label_parsed, strip.position = "left") +
  coord_cartesian(xlim = c(0.01, 1)) +
  expand_limits(y = 0) + 
  labs(x = expression(B[ESH]/B[init]))
ggsave("Figures/meta/Binit.png", g, height = 5, width = 7)

g <- ggplot(SH_order, aes(`B[ESH]/B[init]`, value)) + 
  geom_smooth() + 
  geom_text(aes(label = St), size = 3) + 
  #geom_point() +
  #ggrepel::geom_text_repel(aes(label = St)) + 
  #geom_hline(yintercept = 0.001) + 
  theme_bw() + 
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA, colour = NA), 
        axis.title.y = element_blank(),
        strip.placement = "outside") + 
  facet_wrap(vars(variable), scales = "free_y", labeller = label_parsed, strip.position = "left") +
  coord_cartesian(xlim = c(0.01, 1)) +
  expand_limits(y = 0) + 
  labs(x = expression(B[ESH]/B[init]))
ggsave("Figures/meta/Binit_smooth.png", g, height = 5, width = 7)


SH_order_cor <- SH_order %>% group_by(variable) %>%
  summarise(corr = cor(`B[ESH]/B[init]`, value) %>% round(2), 
            p.value = cor.test(`B[ESH]/B[init]`, value)$p.value) %>%
  #mutate(label = ifelse(p.value <= 0.05, paste0(corr, "*"), corr)) %>%
  mutate(label = corr)


g <- ggplot(SH_order, aes(`B[ESH]/B[init]`, value)) + 
  #geom_smooth() + 
  geom_text(aes(label = St, colour = as.factor(tv)), size = 3) + 
  geom_label(data = SH_order_cor, 
            aes(label = label), 
            x = Inf, 
            y = Inf, 
            vjust = "inward", 
            hjust = "inward") + 
  theme_bw() + 
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA, colour = NA), 
        axis.title.y = element_blank(),
        strip.placement = "outside") + 
  facet_wrap(vars(variable), 
             scales = "free_y", 
             labeller = label_parsed,
             strip.position = "left") +
  expand_limits(x = 0, y = 0) + 
  labs(x = expression(B[ESH]/B[init])) +
  scale_colour_manual(values = c("0" = "black", "1" = "red"))
ggsave("Figures/meta/Binit_tv.png", g, height = 5, width = 7)







SH_pairs <- reshape2::acast(SH_order, list("St", "variable"))

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  val <- cor(x, y, use = "complete.obs")
  r <- abs(val)
  txt <- format(c(val, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  #if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * r)
  text(0.5, 0.5, txt, cex = 2)
}

text2 <- function(x, y, ...) {
  text(x, y, labels = rownames(SH_pairs), col = ifelse(SH$tv, "red", "black"),
       xlim = c(0, 1.1 * max(x)), ylim = c(0, 1.1 * max(y)))
  if(cor(x, y, use = "complete.obs") > 0) abline(a = 0, b = 1, lty = 3)
}

png("Figures/meta/Fpairs_tv.png", height = 4, width = 4, units = "in", res = 400)
pairs(SH_pairs[, 1:3], 
      panel = text2,
      labels = colnames(SH_pairs)[1:3] %>% parse(text = .),
      gap = 0,
      lower.panel = panel.cor)
dev.off()

png("Figures/meta/Bpairs_tv.png", height = 6, width = 6, units = "in", res = 400)
pairs(SH_pairs[, -c(1:3)], 
      panel = text2,
      labels = colnames(SH_pairs)[-c(1:3)] %>% parse(text = .),
      gap = 0,
      lower.panel = panel.cor)
dev.off()


