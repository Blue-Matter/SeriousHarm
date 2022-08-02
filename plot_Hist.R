
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
    mutate(YESH = YESH, dep = dep) %>%
    select(Year, SB, Stock, YESH, dep)
}) %>% bind_rows()

lev_stock <- local({
  x <- group_by(SBall, Stock) %>% summarise(YESH = unique(YESH), dep = unique(dep))
  Stock <- x$Stock[order(x$dep)]
  
  x$Stock <- factor(x$Stock, levels = x$Stock[order(x$dep)])
  x
})

SBall$Stock <- SBall$Stock %>% factor(levels = levels(lev_stock$Stock))


g <- ggplot(SBall, aes(Year, SB)) + 
  geom_line() +
  theme_bw() + 
  geom_vline(data = lev_stock, aes(xintercept = YESH), linetype = 2) + 
  facet_wrap(vars(Stock),
             ncol = 4,
             scales = "free",
             labeller = label_parsed) +
  expand_limits(y = 0) +
  theme(strip.background = element_rect(fill = NA, colour = NA),
        panel.spacing = unit(0, "in"))
ggsave("Figures/meta/SB.png", g, height = 4)



for(i in stocks) {
  
  LRP <- readRDS(file = paste0("LRP/LRP_", i, ".rds"))
  
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
  
  message(i)
  
  g1 <- gghist(LRP, title = TRUE, letter.legend = paste0("(", letters[1:3], ")"))
  
  g2 <- ggSR_LRP(LRP, title = FALSE, panel_labs = c("SB[0.5~Rmax]", "SB[0.9~R/S]"),
                 letter.legend = paste0("(", letters[4:5], ")"))
  
  g3 <- ggSP(LRP, SPB = FALSE, title = FALSE, letter.legend = "(f)")
  
  gright <- cowplot::plot_grid(g2, g3, ncol = 1)
  gout <- cowplot::plot_grid(g1, gright)
  
  ggsave(paste0("Figures/case_studies/", i, ".png"), gout, height = 5, width = 8)
  
}
