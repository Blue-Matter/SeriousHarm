



# Plot all historical figures
source("99-analysis_fn.R")


################# Calculate within-stock correlations of LRPs
stock_order <- readr::read_csv("Tables/stock_order.csv")

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  
  if(!all(is.na(x)) && !all(is.na(y))) {
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
  
}

panel.phase <- function(x, y, yr, pos.txt, ...) {
  if(!all(is.na(x)) && !all(is.na(y))) {
    color_manual <- viridisLite::viridis(length(x) - 1)
    segments(x[2:length(x) - 1], y[2:length(x) - 1], x[2:length(x)], y[2:length(x)], col = color_manual, lwd = 2)
    points(x[c(1, length(x))], y[c(1, length(x))], pch = c(1, 16), cex = 1.5)
    text(x[c(1, length(x))], y[c(1, length(x))], labels = yr, pos = pos.txt)
  }
}


for(j in 1:nrow(stock_order)) {
  i <- stock_order$code[j]
  n <- stock_order$n[j]
  
  LRP <- readRDS(file = paste0("LRP/LRP_", i, ".rds"))
  
  F_LRP <- LRP[[2]] %>% filter(Year <= LRP$Assess$Year) %>% 
    transmute(Year = Year, `F/F[MSY]` = FM/FMSY, `F/F[rep]` = FM/Fmed, `SPR[eq]` = SPR, `SPR[dyn]` = SPRdyn)
  
  png(paste0("Figures/within/FLRP_", n, "_", i, ".png"), width = 6, height = 6, units = "in", res = 400)
  pairs(F_LRP[, -1], 
        labels = colnames(F_LRP[, -1]) %>% parse(text = .),
        gap = 0,
        upper.panel = panel.phase, yr = range(F_LRP$Year), pos.txt = 4,
        lower.panel = panel.cor,
        main = parse(text = stock_order$Stock2[j]))
  dev.off()
  
  B_LRP <- LRP[[2]] %>% filter(Year <= LRP$Assess$Year) %>% 
    transmute(Year = Year, `B/B[MSY]` = SB/SBMSY, `B/B[0~eq]` = SB/SB0, `B/B[0~dyn]` = SB/SB0dyn, 
              `B/B[recover]` = SB/SB_recover, `B/B[SRP]` = SB/SB_SRP, `B/B[0.5~Rmax]` = SB/SB_Rmax, `B/B[0.9~R/S]` = SB/SB_90)
  
  no_LRP <- apply(B_LRP, 2, function(x) all(is.na(x)))
  B_LRP <- B_LRP[, !no_LRP]
  
  png(paste0("Figures/within/BLRP_", n, "_", i, ".png"), width = 6, height = 6, units = "in", res = 400)
  pairs(B_LRP[, -1], 
        labels = colnames(B_LRP[, -1]) %>% parse(text = .),
        gap = 0,
        upper.panel = panel.phase, yr = range(B_LRP$Year), pos.txt = 4,
        lower.panel = panel.cor,
        main = parse(text = stock_order$Stock2[j]))
  dev.off()
}
