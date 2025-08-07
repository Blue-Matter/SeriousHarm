
library(MSEtool)

source("99-analysis_fn.R")


# Plot all historical figures
stocks <- c("ICCAT_SWO", "GoM_RS", "NE_redfish", "sne_yt", "SA_RedPorgy", 
            "SCB_cowcod", "NEA_hom", 
            "WC_darkblotched", "sGSL_cod", "sGSL_herring", "NAFO_plaice", "NAFO_cod", "GoM_cod",
            "ICCAT_BET", "WCVI_herring", #"EBS_cod", 
            "SA_gag")

# Plot SB for all
SBall <- lapply(stocks, function(i) {
  LRP <- readRDS(file = paste0("LRP/LRP_", i, ".rds"))
  YESH <- LRP$Assess$Year
  dep <- LRP$RPts$SB[LRP$RPts$Year == YESH]/LRP$RPts$SB[1]
  LRP$RPts %>% 
    select(Year, SB, Stock, Removals) %>%
    mutate(YESH = YESH, dep = round(dep, 2), code = i)
}) %>% bind_rows()

stock_order <- SBall %>%
  group_by(Stock) %>%
  summarise(dep = mean(dep), code = unique(code), `Year ESH` = unique(YESH)) %>%
  arrange(dep) %>%
  mutate(n = 1:n(),
         Stock2 = paste0("\"(", n, ")\"~", Stock))
readr::write_csv(stock_order, "Tables/stock_order.csv")

# Update stock summary table
stock_sumry <- left_join(
  stock_order %>% select(!Stock),
  readr::read_csv("Tables/stock_summary.csv")
) %>%
  mutate(Stock = paste0("(", n, ") ", Stock)) %>%
  select(Stock, Region, `Time-varying productivity`, `Assessment Time Period`, `Stock-recruit relationship`, `Year ESH`, `ESH criteria`, References)
readr::write_csv(stock_sumry, "Tables/stock_summary_update.csv")

# Make SB figures
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


#g <- make_SB_figure(SBout, stocks = stock_order$Stock) %>%
#  ggpubr::annotate_figure(bottom = "Year")
#ggsave("Figures/meta/SB_with_catch.png", g, height = 8, width = 7)

g <- make_SB_figure(SBout, stocks = stock_order$Stock, lim = "all") %>%
  ggpubr::annotate_figure(bottom = "Year")
ggsave("Figures/meta/SB_with_catch_all.png", g, height = 8, width = 7)

# All case studies
for(j in 1:nrow(stock_order)) {
  
  i <- stock_order$code[j]
  n <- stock_order$n[j]
  
  LRP <- readRDS(file = paste0("LRP/LRP_", i, ".rds"))
  LRP[["RPvars"]][["Stock"]] <- LRP[["RPts"]]$Stock <- stock_order %>% filter(code == i) %>% pull(Stock2)
  
  message(i)
  message("Year range: ", range(LRP$RPts$Year) %>% paste0(collapse = " - "))
  message("Year of ESH: ", LRP$Assess$Year)
  message("B_ESH/Binit = ", round(LRP$RPts$SB[LRP$RPts$Year == LRP$Assess$Year]/LRP$RPts$SB[1], 2))
  message("Time-varying phi0: ", length(unique(LRP$RPts$phi0)) > 1)
    
  g1 <- plot_Hist(LRP)
  ggsave(paste0("Figures/Hist/ref_", n, "_", i, ".png"), g1, height = 5, width = 6)
  
  if(!is.na(LRP$RPvars$Yr_SP)) message("Low SP for ", i)
  g2 <- plot_SP(LRP)
  ggsave(paste0("Figures/Hist/SP_", n, "_", i, ".png"), g2, height = 4, width = 4)
  
  g3 <- plot_SR_LRP(LRP, year = FALSE)
  ggsave(paste0("Figures/Hist/SRR_", n, "_", i, ".png"), g3, height = 4, width = 4)
  
  g4 <- plot_ts(LRP)
  ggsave(paste0("Figures/Hist/ts_", n, "_", i, ".png"), g4, height = 4, width = 8)
  
}

# Case study figures
#case_studies <- c("sne_yt", "sGSL_herring", "WCVI_herring", "GoM_RS")
case_studies <- stocks

for (j in 1:nrow(stock_order)) {
  
  i <- stock_order$code[j]
  n <- stock_order$n[j]
  
  LRP <- readRDS(file = paste0("LRP/LRP_", i, ".rds"))
  LRP[["RPvars"]][["Stock"]] <- LRP[["RPts"]]$Stock <- stock_order %>% filter(code == i) %>% pull(Stock2)
  
  message(i)
  
  g1 <- gghist(LRP, title = TRUE, letter.legend = paste0("(", letters[1:3], ")"))
  
  g2 <- ggSR_LRP(LRP, title = FALSE, panel_labs = c("SB[0.5~Rmax]", "SB[0.9~R/S]"),
                 letter.legend = paste0("(", letters[4:5], ")"))
  
  g3 <- ggSP(LRP, SPB = TRUE, title = FALSE, letter.legend = "(f)")
  
  gright <- cowplot::plot_grid(g2, g3, ncol = 1)
  gout <- cowplot::plot_grid(g1, gright)
  
  ggsave(paste0("Figures/case_studies/", n, "_", i, ".png"), gout, height = 5, width = 8)
}


# Calculate reference points
ref_table <- lapply(1:length(stocks), function(i) {
  
  ref <- readRDS(file.path("LRP", paste0("LRP_", stocks[i], ".rds")))
  
  ref$Assess <- round(ref$Assess, 2)
  
  out <- tibble::tibble(
    Stock = ref$RPvars$Stock,
    `F/F[MSY]` = ref$Assess$FMSY,
    `SPR[eq]` = ref$Assess$SPR,
    `SPR[dyn]` = ref$Assess$SPRdyn,
    `SB/SB[MSY]` = ref$Assess$SBMSY,
    `SB/SB[0,eq]` = ref$Assess$SB0,
    `SB/SB[0,dyn]` = ref$Assess$SB0dyn,
    `SB/SB[recover]` = ref$Assess$SB_recover,
    `SB/SB[lowSP]` = ref$Assess$SB_SP,
    `SB/SB[0.5Rmax]` = ref$Assess$SB_Rmax,
    `SB/SB[0.9R/S]` = ref$Assess$SB_90
  ) %>%
    left_join(stock_order)
  
  if (is.na(out$`SB/SB[lowSP]`)) {
    out$`SB/SB[lowSP]` <- "-"
  } else {
    out$`SB/SB[lowSP]` <- paste0(out$`SB/SB[lowSP]`, " (", ref$RPvars$Yr_SP, ")")
  }
  
  if (is.na(out$`SB/SB[recover]`)) {
    out$`SB/SB[recover]` <- "-"
  } else {
    out$`SB/SB[recover]` <- paste0(out$`SB/SB[recover]`, " (", ref$RPvars$Yr_Brecover, ")")
  }
  
  # For plaice
  #if (is.infinite(out$`F/F[MSY]`)) out$`F/F[MSY]` <- "-"
  #if (is.infinite(out$`SB/SB[MSY]`)) out$`SB/SB[MSY]` <- "-"
  #if (is.infinite(out$`SB/SB[0,eq]`)) out$`SB/SB[0,eq]` <- "-"
  
  return(out)
}) %>%
  bind_rows()

ref_F <- ref_table[order(ref_table$n), ] %>%
  mutate(Stock = paste0("(", n, ") ", gsub("~", " ", Stock))) %>%
  select(
    Stock,
    dep,
    `Year ESH`,
    `F/F[MSY]`,
    `SPR[eq]`,
    `SPR[dyn]`,
    n
  )
readr::write_csv(ref_F, "Tables/LRP_F.csv")


ref_B <- ref_table[order(ref_table$n), ] %>%
  mutate(Stock = paste0("(", n, ") ", gsub("~", " ", Stock))) %>%
  select(
    Stock,
    dep,
    `Year ESH`,
    `SB/SB[MSY]`,
    `SB/SB[0,eq]`,
    `SB/SB[0,dyn]`,
    `SB/SB[recover]`,
    `SB/SB[lowSP]`,
    `SB/SB[0.5Rmax]`,
    `SB/SB[0.9R/S]`,
    n
  )
readr::write_csv(ref_B, "Tables/LRP_B.csv")


ref_table_num <- lapply(1:length(stocks), function(i) {
  
  ref <- readRDS(file.path("LRP", paste0("LRP_", stocks[i], ".rds")))
  
  ref$Assess <- round(ref$Assess, 2)
  
  out <- tibble::tibble(
    Stock = ref$RPvars$Stock,
    `F/F[MSY]` = ref$Assess$FMSY,
    `SPR[eq]` = ref$Assess$SPR,
    `SPR[dyn]` = ref$Assess$SPRdyn,
    `SB/SB[MSY]` = ref$Assess$SBMSY,
    `SB/SB[0,eq]` = ref$Assess$SB0,
    `SB/SB[0,dyn]` = ref$Assess$SB0dyn,
    `SB/SB[recover]` = ref$Assess$SB_recover,
    `SB/SB[lowSP]` = ref$Assess$SB_SP,
    `SB/SB[0.5Rmax]` = ref$Assess$SB_Rmax,
    `SB/SB[0.9R/S]` = ref$Assess$SB_90
  ) %>%
    left_join(stock_order)
  
  return(out)
}) %>%
  bind_rows()

readr::write_csv(ref_table_num, "Tables/LRP_numeric.csv")
