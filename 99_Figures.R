

############# Phase plot of true vs estimated BMSY
library(RPC)
library(MSEtool)
OM <- DFO_Pacific_Cod_2020 %>% SubCpars(1:3)
Hist <- runMSE(OM, Hist = TRUE)
res <- SCA(Data = Hist@Data, AddInd = 1:5)

# SCA_8040
SCA_pcod <- function(x, Data, ...) {
  if(max(Data@Year) > Data@LHYear) {
    Data@AddInd[, 1:4, match(Data@LHYear + 1, Data@Year):length(Data@Year)] <- NA_real_
  }
  SCA(x, Data, AddInd = 1:5)
}
class(SCA_pcod) <- "Assess"
SCA_8040MSY <- make_MP(SCA_pcod, HCR_ramp, OCP_type = "SSB_SSBMSY", LOCP = 0.4, TOCP = 0.8, AddInd = 1:5, diagnostic = "full")
MSE <- Project(Hist, MPs = c("NFref", "SCA_8040MSY"))
saveRDS(MSE, file = "99_HCR_Figure.rds")


MSE <- readRDS("99_HCR_Figure.rds")

sim <- 2
BOM <- data.frame(Year = MSE@OM$CurrentYr[1] + 1:MSE@proyears, OM = MSE@SB_SBMSY[sim, 2, ])
Best <- lapply(MSE@PPD[[2]]@Misc[[sim]]$Assessment_report, function(x) {
  data.frame(Year = names(x@SSB_SSBMSY)[length(x@SSB_SSBMSY)] %>% as.numeric(),
             Est = x@SSB_SSBMSY[length(x@SSB_SSBMSY)])
}) %>% bind_rows()

g1 <- left_join(BOM, Best, by = "Year") %>% dplyr::filter(!is.na(Est)) %>%
  ggplot(aes(OM, Est, colour = Year, label = Year)) + 
  geom_hline(yintercept = c(0.4, 0.8), linetype = 3) + 
  geom_vline(xintercept = c(0.4, 0.8), linetype = 3) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  geom_point() + 
  geom_path() + ggrepel::geom_text_repel() +
  theme_bw() + coord_cartesian(xlim = c(0, 3), ylim = c(0, 3)) + 
  scale_colour_viridis_c() +
  labs(x = expression(True~B/B[MSY]~"(for evaluating performance)"), y = expression(Estimated~B/B[MSY]~"(for HCR)"))
ggsave("Estimated_versus_true_status.png", g1, height = 4, width = 5.5)


FOM <- data.frame(Year = MSE@OM$CurrentYr[1] + 1:MSE@proyears, OM = MSE@F_FMSY[sim, 2, ])
Fest <- lapply(MSE@PPD[[2]]@Misc[[1]]$Assessment_report, function(x) {
  Catch <- SAMtool::HCR_ramp(x, OCP_type = "SSB_SSBMSY", LOCP = 0.4, TOCP = 0.8)@TAC
  
  
  data.frame(Year = names(x@F_FMSY)[length(x@F_FMSY)] %>% as.numeric(),
             Est = x@SSB_SSBMSY[length(x@SSB_SSBMSY)])
}) %>% bind_rows()

relB <- seq(0, 4, 0.2)

g2 <- BOM %>% dplyr::left_join(FOM, by = "Year") %>% 
  dplyr::filter(Year %in% Best$Year) %>%
  mutate(F_HCR = SAMtool::HCRlin(OM.x, 0.4, 0.8, 0, 1)) %>%
  ggplot(aes(OM.x, OM.y, colour = Year, label = Year)) + 
  geom_line(data = data.frame(OM.x = relB,
                              OM.y = ifelse(relB < 0.4, 0, 
                                            ifelse(relB > 1, 1,
                                            SAMtool::HCRlin(relB, 0.4, 0.8, 0, 1)))),
            aes(x = OM.x, y = OM.y), inherit.aes = FALSE) + 
  geom_point() + 
  geom_path(aes(colour = Year)) + 
  ggrepel::geom_text_repel() +
  theme_bw() + coord_cartesian(xlim = c(0, 3), ylim = c(0, 1.75)) + 
  scale_colour_viridis_c() +
  labs(x = expression(True~B/B[MSY]), y = expression(True~F/F[MSY]))
ggsave("Implemented_vs_true_HCR.png", g2, height = 4, width = 5.5)


g3 <- Best %>% mutate(F_HCR = SAMtool::HCRlin(Est, 0.4, 0.8, 0, 1)) %>%
  ggplot(aes(Est, F_HCR, colour = Year, label = Year)) + 
  geom_line(data = data.frame(OM.x = relB,
                              OM.y = ifelse(relB < 0.4, 0, 
                                            ifelse(relB > 1, 1,
                                                   SAMtool::HCRlin(relB, 0.4, 0.8, 0, 1)))),
            aes(x = OM.x, y = OM.y), inherit.aes = FALSE) + 
  geom_point() + 
  ggrepel::geom_text_repel() +
  theme_bw() + coord_cartesian(xlim = c(0, 3), ylim = c(0, 1.75)) + 
  scale_colour_viridis_c() +
  labs(x = expression(Estimated~B/B[MSY]), y = expression(Prescribed~F/F[MSY]))
ggsave("Implemented_HCR.png", g3, height = 4, width = 5.5)

ggpubr::ggarrange(g1, g3, g2, ncol = 1, common.legend = TRUE, legend = "right", align = "hv")
ggsave("HCR.png", height = 7, width = 4)


#################### Single plot that compares P(SH) vs. F/M for 2 stocks

Pcod_PM <- local({
  
  source("00_performance_measures.R")
  
  # Show historical probability of serious harm
  MSE <- readRDS("Pcod/Pcod_M.rds")
  
  type <- c("HistSSB", "iSCAM_SSB0", "SSBMSY", "depletion", "50%Rmax", "90%RS")
  frac <- c(1, 0.3, 0.4, 0.3, 1, 1)
  HistSSB_y <- 2000
  
  PM_hist <- Map(function(x, y) SHPM(MSE, type = x, frac = y, HistSSB_y = HistSSB_y, dep_type = "Dynamic", Hist = TRUE), x = type, y = frac)
  
  annual_PM_hist <- lapply(PM_hist, function(x) data.frame(Year = seq(MSE@OM$CurrentYr[1] - MSE@nyears + 1, MSE@OM$CurrentYr[1]),
                                                           value = x@Prob, Metric = x@Caption)) %>%
    dplyr::bind_rows() %>% mutate(Stock = "Pacific cod")
  
  # Probability of being above 5 definitions of serious harm from projections
  PM <- Map(function(x, y) SHPM(MSE, type = x, frac = y, HistSSB_y = HistSSB_y, dep_type = "Dynamic"), x = type, y = frac)
  
  annual_PM <- lapply(PM, function(x) {
    apply(x@Stat >= x@Ref, 2:3, mean) %>% structure(dimnames = list(FM = relF, Year = 2020 + 1:MSE@proyears)) %>% reshape2::melt() %>% 
      dplyr::mutate(Metric = x@Caption, Name = x@Name, FMlabel = paste0("F/M==", FM))
  }) %>% dplyr::bind_rows() %>% mutate(Stock = "Pacific cod")
  
  list(annual_PM_hist, annual_PM)
})

Pherring_PM <- local({
  
  source("00_performance_measures.R")
  
  MSE <- readRDS("Pherring/Pherring_M.rds")
  
  # Show historical probability of serious harm
  type <- c("HistSSB", "iSCAM_SSB0", "SSBMSY", "depletion", "50%Rmax", "90%RS")
  frac <- c(1, 0.3, 0.4, 0.3, 1, 1)
  HistSSB_y <- 2010
  
  PM_hist <- Map(function(x, y) SHPM(MSE, type = x, frac = y, HistSSB_y = HistSSB_y, dep_type = "Dynamic", Hist = TRUE), x = type, y = frac)
  
  annual_PM_hist <- lapply(PM_hist, function(x) data.frame(Year = seq(MSE@OM$CurrentYr[1] - MSE@nyears + 1, MSE@OM$CurrentYr[1]),
                                                           value = x@Prob, Metric = x@Caption)) %>%
    dplyr::bind_rows() %>% mutate(Stock = "WCVI Pacific herring")
  
  # Probability of being above 5 definitions of serious harm from projections
  PM <- Map(function(x, y) SHPM(MSE, type = x, frac = y, HistSSB_y = HistSSB_y, dep_type = "Dynamic"), x = type, y = frac)
  
  annual_PM <- lapply(PM, function(x) {
    apply(x@Stat >= x@Ref, 2:3, mean) %>% 
      structure(dimnames = list(FM = relF, Year = MSE@OM$CurrentYr[1] + 1:MSE@proyears)) %>% reshape2::melt() %>%
      dplyr::mutate(Metric = x@Caption, Name = x@Name, FMlabel = paste0("F/M==", FM))
  }) %>% dplyr::bind_rows() %>% mutate(Stock = "WCVI Pacific herring")
  
  list(annual_PM_hist, annual_PM)
})


# Show historical probability of serious harm
PM_hist <- rbind(Pherring_PM[[1]], Pcod_PM[[1]])
ggplot(PM_hist, aes(Year, value, colour = Metric)) + facet_grid(Stock ~ ., scales = "free_x") + geom_point() + geom_line() + theme_bw() + 
  scale_colour_discrete(labels = scales::label_parse()) +
  labs(y = "Probability above serious harm", colour = "Serious harm candidates")
ggsave("Serious_harm_Hist.png", height = 5, width = 7)

PM <- rbind(Pherring_PM[[2]], Pcod_PM[[2]]) %>% group_by(Stock) %>% dplyr::filter(Year == max(Year))

# Figure
ggplot(PM, aes(FM, value, colour = Metric)) + facet_grid(Stock ~ .) + geom_point() + geom_line() + theme_bw() + 
  scale_colour_discrete(labels = scales::label_parse()) +
  labs(x = expression(F/M), y = "Probability above serious harm after 50 years", colour = "Serious harm candidates")
ggsave("Serious_harm_FM.png", height = 4, width = 6)

# Table (FM = 0)
PM_table <- PM %>% dplyr::filter(FM == 0) %>% mutate(value = round(value, 2)) %>% 
  reshape2::dcast(list("Stock", "Name"), variable.name = "value") %>%
  structure(names = c("Stock", "30% Dynamic SSB0", "30% iSCAM SSB0", "SSB_50%Rmax", "SSB_90%R/S", "SSB_2000", "SSB_2010", "40% SSBMSY"))
write.csv(PM_table, file = "Serious_harm_FM.csv")


# FM ratio such that P(B > SH) = 0.95
LRP_levels <- c("30% Dynamic SSB0", "30% iSCAM SSB0", "SSB_50%Rmax", "SSB_90%R/S", "SSB_2000", "SSB_2010", "40% SSBMSY")

Pcod_FP95 <- local({
  type <- c("HistSSB", "iSCAM_SSB0", "SSBMSY", "depletion", "50%Rmax", "90%RS")
  FSH <- readRDS("Pcod/Pcod_FP95_M.rds") %>% sapply(function(x) ifelse(x$FSH <= 0.01, 0, x$FSH))
  data.frame(FP95 = FSH[c(4, 2, 5, 6, 1, NA, 3)], Stock = "Pacific cod", LRP = factor(LRP, levels = LRP))
})

Pherring_FP95 <- local({
  type <- c("HistSSB", "iSCAM_SSB0", "SSBMSY", "depletion", "50%Rmax", "90%RS")
  FSH <- readRDS("Pherring/Pherring_FP95_M.rds") %>% sapply(function(x) ifelse(x$FSH <= 0.01, 0, x$FSH))
  data.frame(FP95 = FSH[c(4, 2, 5, 6, NA, 1, 3)], Stock = "WCVI Pacific herring", LRP = factor(LRP, levels = LRP))
})

rbind(Pcod_FP95, Pherring_FP95) %>% reshape2::dcast(list("Stock", "LRP"), value.var = "FP95") %>%
  write.csv(file = "FP95.csv")
