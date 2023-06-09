
library(MSEtool)

source("analysis_fn.R")


# Plot all historical figures
stocks <- c("ICCAT_SWO", "GoM_RS", "NE_redfish", "sne_yt", "SA_RedPorgy", "SCB_cowcod", "NEA_hom", 
            "WC_darkblotched", "sGSL_cod", "sGSL_herring", "NAFO_plaice", "NAFO_cod", "GoM_cod",
            "ICCAT_BET", "WCVI_herring", "EBS_cod", "SA_gag")

# Plot SB for all
SBall <- lapply(stocks, function(i) {
  LRP <- readRDS(file = paste0("LRP/LRP_", i, ".rds"))
  YESH <- LRP$Assess$Year
  dep <- LRP$RPts$SB[LRP$RPts$Year == YESH]/LRP$RPts$SB[1]
  LRP$RPts %>% 
    select(Year, SB, Stock, Removals) %>%
    mutate(YESH = YESH, dep = dep, code = i)
}) %>% bind_rows()

stock_order <- SBall %>%
  group_by(Stock) %>%
  summarise(dep = mean(dep), code = unique(code)) %>%
  arrange(dep) %>%
  mutate(n = 1:n(),
         Stock2 = paste0("\"(", n, ")\"~", Stock)) %>%
  select(-dep)

SBout <- left_join(SBall, stock_order, by = c("Stock", "code"))

make_SB_figure <- function(SBout, stocks, lim = c("YESH", "all")) {
  lim <- match.arg(lim)
  
  g <- lapply(stocks, function(sp) {
    dframe <- SBout %>%
      filter(Stock == sp) %>%
      select(Year, SB, Removals, Stock2, YESH) %>%
      rename(Catch = Removals) %>%
      reshape2::melt(id.vars = c("Year", "Stock2", "YESH")) %>%
      mutate(SB = ifelse(variable == "SB", value, NA_real_),
             Catch = ifelse(variable == "SB", NA_real_, value))
    
    dframe_yesh <- filter(dframe, Year <= YESH)
    dframe_future <- filter(dframe, Year > YESH)
    g_out <- dframe_yesh %>%
      ggplot(aes(x = Year)) + 
      geom_line(aes(y = SB)) +
      geom_col(colour = "black", linewidth = 0.1, fill = "grey40", width = 1, aes(y = Catch)) + 
      geom_point(data = dframe %>% filter(Year == YESH, variable == "SB"), aes(y = SB)) +
      facet_grid(vars(variable), 
                 vars(Stock2), 
                 scales = "free_y", 
                 switch = "y",
                 labeller = label_parsed) +
      expand_limits(y = 0) +
      #labs(x = "", y = "") +
      theme_bw() + 
      theme(strip.background = element_rect(fill = NA, colour = NA), 
            axis.title = element_blank(),
            panel.spacing = unit(0, "in"),
            strip.text.y = element_blank(),
            axis.text.y = element_blank(),
            #axis.text.x = element_text(angle = 45),
            axis.ticks.length.y = unit(0, "in"))
    
    if (lim == "all") {
      g_out <- g_out +
        geom_line(data = dframe_future, aes(y = SB), linetype = 4) +
        geom_col(data = dframe_future, aes(y = Catch), width = 1, linewidth = 0.1, colour = "black", fill = "white")
    }
    g_out
  })
  cowplot::plot_grid(plotlist = g, ncol = 4, align = "hv")
}


g <- make_SB_figure(SBout, stocks = stock_order$Stock[1:2]) %>%
  ggpubr::annotate_figure(bottom = "Year")
ggsave("Figures/meta/SB_with_catch.png", g, height = 8, width = 7)

g <- make_SB_figure(SBout, stocks = stock_order$Stock[1:2], lim = "all") %>%
  ggpubr::annotate_figure(bottom = "Year")
ggsave("Figures/meta/SB_with_catch_all.png", g, height = 8, width = 7)

# All case studies
for(i in stocks) {
  
  LRP <- readRDS(file = paste0("LRP/LRP_", i, ".rds"))
  LRP[["RPvars"]][["Stock"]] <- LRP[["RPts"]]$Stock <- stock_order %>% filter(code == i) %>% pull(Stock2)
  
  message(i)
  message("Year range: ", range(LRP$RPts$Year) %>% paste0(collapse = " - "))
  message("Year of ESH: ", LRP$Assess$Year)
  message("B_ESH/Binit = ", round(LRP$RPts$SB[LRP$RPts$Year == LRP$Assess$Year]/LRP$RPts$SB[1], 2))
  message("Time-varying phi0: ", length(unique(LRP$RPts$phi0)) > 1)
    
  g1 <- plot_Hist(LRP)
  ggsave(paste0("Figures/Hist/", i, ".png"), g1, height = 6, width = 4)
  
  if(!is.na(LRP$RPvars$Yr_SP)) message("Low SP for ", i)
  g2 <- plot_SP(LRP)
  ggsave(paste0("Figures/SP/SP_", i, ".png"), g2, height = 6, width = 4)
  
  g3 <- plot_SR_LRP(LRP, year = FALSE)
  ggsave(paste0("Figures/SRR/SRR_", i, ".png"), g3, height = 4, width = 4)
  
  g4 <- plot_ts(LRP)
  ggsave(paste0("Figures/ts/ts_", i, ".png"), g4, height = 4, width = 8)
  
}

# Case study figures
#case_studies <- c("sne_yt", "sGSL_herring", "WCVI_herring", "GoM_RS")
case_studies <- stocks

for(i in case_studies) {
  
  LRP <- readRDS(file = paste0("LRP/LRP_", i, ".rds"))
  LRP[["RPvars"]][["Stock"]] <- LRP[["RPts"]]$Stock <- stock_order %>% filter(code == i) %>% pull(Stock2)
  
  message(i)
  
  g1 <- gghist(LRP, title = TRUE, letter.legend = paste0("(", letters[1:3], ")"))
  
  g2 <- ggSR_LRP(LRP, title = FALSE, panel_labs = c("SB[0.5~Rmax]", "SB[0.9~R/S]"),
                 letter.legend = paste0("(", letters[4:5], ")"))
  
  g3 <- ggSP(LRP, SPB = TRUE, title = FALSE, letter.legend = "(f)")
  
  gright <- cowplot::plot_grid(g2, g3, ncol = 1)
  gout <- cowplot::plot_grid(g1, gright)
  
  ggsave(paste0("Figures/case_studies/", i, ".png"), gout, height = 5, width = 8)
}
