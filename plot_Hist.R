
library(MSEtool)
source("analysis_fn.R")


# Plot all historical figures
stocks <- c("ICCAT_SWO", "GoM_RS", "NE_redfish", "GoM_haddock", "SA_RedPorgy", "SCB_cowcod", "NEA_hom", 
            "WC_darkblotched", "WC_pop", "sGSL_cod", "sGSL_herring", "NAFO_plaice", "NAFO_cod", "GoM_cod",
            "ICCAT_BET", "WCVI_herring", "EBS_cod")


for(i in stocks) {
  
  LRP <- readRDS(file = paste0("LRP/LRP_", i, ".rds"))
  
  message(i)
  message("Year range: ", range(LRP$RPts$Year) %>% paste0(collapse = " - "))
  message("Year of ESH: ", LRP$Assess$Year)
  message("B_ESH/Binit = ", round(LRP$RPts$SB[LRP$RPts$Year == LRP$Assess$Year]/LRP$RPts$SB[1], 2))
  message("Time-varying phi0: ", length(unique(LRP$RPts$phi0)) > 1)
  
  g <- plot_Hist(LRP)
  ggsave(paste0("Figures/Hist/", i, ".png"), g, height = 6, width = 4)
  
  g <- plot_SP(LRP)
  ggsave(paste0("Figures/SP/SP_", i, ".png"), g, height = 6, width = 4)
  
  g <- plot_SR_LRP(LRP, year = TRUE)
  ggsave(paste0("Figures/SRR/SP_", i, ".png"), g, height = 4, width = 4)
  
  g <- plot_ts(LRP)
  ggsave(paste0("Figures/ts/ts_", i, ".png"), g, height = 4, width = 8)
}

