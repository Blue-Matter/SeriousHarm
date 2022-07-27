
library(MSEtool)
source("analysis_fn.R")


# Plot all historical figures
stocks <- c("ICCAT_SWO", "GoM_RS", "NE_redfish", "sne_yt", "SA_RedPorgy", "SCB_cowcod", "NEA_hom", 
            "WC_darkblotched", "sGSL_cod", "sGSL_herring", "NAFO_plaice", "NAFO_cod", "GoM_cod",
            "ICCAT_BET", "WCVI_herring", "EBS_cod", "SA_gag")



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

