


library(tidyverse)

SH <- readxl::read_excel('Report/Seriousharm_estimates_import.xlsx', sheet = "Seriousharm_estimates") %>% 
  select(!`F/F[rep]` & !`SB/SB[SRP]`) %>%
  filter(Stock != "US Pacific ocean perch", Stock != "GM haddock")

SH <- SH[order(SH$`B[ESH]/B[init]`), ] %>% 
  #mutate(St = paste0("(", 1:nrow(SH), ")")) %>% 
  mutate(St = 1:nrow(SH))



###### For ggplot
SH_order <- SH %>% 
  reshape2::melt(id.vars = c("Stock", "St", "Year_ESH", "B[ESH]/B[init]", "tv")) %>%
  mutate(value_out = ifelse(value > 3, 3, value),
         St = ifelse(value > 3, paste0(St, "*"), St)) %>%
  filter(!is.na(value))
var_lev <- unique(SH_order$variable)
levels(SH_order$variable) <- var_lev[c(1:6, 9, 8, 7, 10)]

SH_order_cor <- SH_order %>% group_by(variable) %>%
  summarise(corr = cor(`B[ESH]/B[init]`, value) %>% round(2), 
            p.value = cor.test(`B[ESH]/B[init]`, value)$p.value) %>%
  #mutate(label = ifelse(p.value <= 0.05, paste0(corr, "*"), corr)) %>%
  mutate(label = corr)

####### For pairs plot
SH_pairs <- SH %>% 
  reshape2::melt(id.vars = c("Stock", "St", "Year_ESH", "B[ESH]/B[init]", "tv")) %>%
  mutate(value_out = ifelse(value > 3, 3, value),
         variable = factor(variable, levels = var_lev[c(1:6, 9, 8, 7, 10)])) %>%
  reshape2::acast(list("St", "variable"), value.var = "value_out")

# Shade indicates region where serious harm is triggered
shade_SH <- data.frame(variable = colnames(SH_pairs),
                       low = c(1, rep(-Inf, ncol(SH_pairs) - 1)),
                       high = c(Inf, 0.4, 0.4, 0.4, 0.2, 0.2, 1, 1, 1, 1))

shade_gg <- lapply(0:1, function(x) mutate(shade_SH, x = x)) %>% 
  bind_rows() %>%
  mutate(variable = factor(variable, levels = var_lev[c(1:6, 9, 8, 7, 10)]))

SH_pairs_plot <- rbind(SH_pairs, shade_SH$low, shade_SH$high)


g <- ggplot(SH_order, aes(`B[ESH]/B[init]`, value_out)) + 
  geom_ribbon(data = shade_gg, aes(x = x, ymin = low, ymax = high), inherit.aes = FALSE, fill = "#CD5C5C") + 
  geom_text(aes(label = St), size = 3) + 
  geom_text(data = SH_order_cor, 
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
             ncol = 3,
             strip.position = "left") +
  coord_cartesian(xlim = c(0.01, 1)) +
  expand_limits(y = 0) + 
  labs(x = expression(Delta~SB))
ggsave("Figures/meta/Binit.png", g, height = 5, width = 6)

g <- ggplot(SH_order, aes(`B[ESH]/B[init]`, value_out)) + 
  geom_ribbon(data = shade_gg, aes(x = x, ymin = low, ymax = high), inherit.aes = FALSE, fill = "#CD5C5C") + 
  geom_smooth() + 
  geom_text(aes(label = St), size = 3) + 
  theme_bw() + 
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA, colour = NA), 
        axis.title.y = element_blank(),
        strip.placement = "outside") + 
  facet_wrap(vars(variable), 
             scales = "free_y", 
             labeller = label_parsed, 
             ncol = 3,
             strip.position = "left") +
  coord_cartesian(xlim = c(0.01, 1)) +
  expand_limits(y = 0) + 
  labs(x = expression(Delta~SB))
ggsave("Figures/meta/Binit_smooth.png", g, height = 5, width = 6)


g <- ggplot(SH_order, aes(`B[ESH]/B[init]`, value_out)) + 
  geom_ribbon(data = shade_gg, aes(x = x, ymin = low, ymax = high), inherit.aes = FALSE, fill = "#CD5C5C") + 
  geom_text(data = SH_order %>% filter(tv == 0), 
            aes(label = St),
            fontface = 1,
            size = 2.5) +
  geom_text(data = SH_order %>% filter(tv == 1), 
            aes(label = St),
            fontface = 4,
            size = 2.5) + 
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
             ncol = 3,
             strip.position = "left") +
  expand_limits(x = 0, y = 0) + 
  labs(x = expression(Delta~SB))
ggsave("Figures/meta/Binit_tv.png", g, height = 5, width = 6)






# Color ramp
graphics::colorRamp()

#x <- seq(-1, 1, 0.01)
#rho_col <- grDevices::colorRampPalette(c("blue", "grey90", "red"))(length(x))
#plot(x, col = rho_col)

color_ramp <- data.frame(x = seq(-1, 1, by = 0.01)) %>%
  mutate(rho_col = grDevices::colorRampPalette(c("turquoise", "grey90", "#DC143C"))(length(x)))


# Panel function that reports correlation and background color
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor = 1, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  val <- cor(x[1:nrow(SH_pairs)], y[1:nrow(SH_pairs)], use = "complete.obs")
  r <- abs(val)
  txt <- format(c(val, 0.123456789) %>% round(2), digits = digits)[1]
  txt <- paste0(prefix, txt)
  #if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * r)
  col_bg <- color_ramp %>%
    filter(round(x, 2) == as.numeric(txt)) %>%
    pull(rho_col)
  
  rect(0, 0, 1, 1, col = col_bg)
  text(0.5, 0.5, txt, cex = cex.cor)
}

panel.cor.F <- panel.cor.B <- panel.cor
formals(panel.cor.F)$cex.cor <- 2



text_SH <- function(x, y, tv = FALSE, cex = 0.5, ...) {
  #browser()
  xx <- x[1:nrow(SH_pairs)]
  yy <- y[1:nrow(SH_pairs)]
  
  x_range <- range(c(0, xx), na.rm = TRUE)
  y_range <- range(c(0, yy), na.rm = TRUE)
  
  # Red zone
  # Colors from https://www.rapidtables.com/web/color/html-color-codes.html
  green <- "#228B22"    #forestgreen
  yellow <- "#F0E68C"   #khaki
  red <- "#CD5C5C"      #indianred
    
  polygon(x = c(x[18], x[19], x[19], x[18]) %>% pmin(max(x_range)) %>% pmax(min(x_range)), 
          y = c(-1e8, -1e8, 1e8, 1e8) %>% pmin(max(y_range)) %>% pmax(min(y_range)), 
          col = yellow,
          border = NA)
  polygon(x = c(-1e8, 1e8, 1e8, -1e8) %>% pmin(max(x_range)) %>% pmax(min(x_range)), 
          y = c(y[18], y[18], y[19], y[19]) %>% pmin(max(y_range)) %>% pmax(min(y_range)), 
          col = yellow, border = NA)
  polygon(x = c(x[18], x[19], x[19], x[18]) %>% pmin(max(x_range)) %>% pmax(min(x_range)), 
          y = c(y[18], y[18], y[19], y[19]) %>% pmin(max(y_range)) %>% pmax(min(y_range)), 
          col = red,
          border = NA)
  
  if(tv) {
    text(c(0, xx), c(0, yy), 
         labels = c(NA, rownames(SH_pairs)),
         font = ifelse(SH$tv, 4, 1), 
         cex = cex)
  } else {
    text(c(0, xx), c(0, yy), 
         labels = c(NA, rownames(SH_pairs)), 
         cex = cex)
  }
  
  if(cor(xx, yy, use = "complete.obs") > 0) abline(a = 0, b = 1, lty = 3)
}


png("Figures/meta/Fpairs_tv.png", height = 4, width = 4, units = "in", res = 400)
pairs(SH_pairs_plot[, 1:3], 
      panel = text_SH,
      cex = 1,
      tv = TRUE,
      labels = colnames(SH_pairs)[1:3] %>% parse(text = .),
      gap = 0,
      lower.panel = panel.cor.F)
dev.off()

png("Figures/meta/Bpairs_tv.png", height = 6, width = 6, units = "in", res = 400)
pairs(SH_pairs_plot[, -c(1:3)], 
      panel = text_SH,
      cex = 0.75,
      tv = TRUE,
      labels = colnames(SH_pairs)[-c(1:3)] %>% parse(text = .),
      gap = 0,
      lower.panel = panel.cor.B)
dev.off()


