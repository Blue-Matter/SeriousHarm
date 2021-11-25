
library(MSEtool)

nsim <- 3
OM <- RPC::DFO_Pacific_Cod_2020 %>% SubCpars(1:nsim)
OM@interval <- 200
Hist <- runMSE(OM, Hist = TRUE, parallel = TRUE)


##### First step - fixed F MPs to identify suitability of serious harm metrics to fishing

source("00_performance_measures.R")

MPs <- paste0(relF, "_FMSY")
MSE <- Project(Hist, MPs = MPs, extended = TRUE)
saveRDS(MSE, file = "Pcod/Pcod_FMSY.rds")

MSE <- readRDS("Pcod/Pcod_FMSY.rds")

##### First step - fixed F MPs to identify suitability of serious harm metrics to fishing
MPs <- paste0(relF, "_M")
MSE <- Project(Hist, MPs = MPs, extended = TRUE)
saveRDS(MSE, file = "Pcod/Pcod_M.rds")

MSE <- readRDS("Pcod/Pcod_FMSY.rds")


##### Calculate F such that P(SH) = 0.95
source("00_Fopt.R")

type <- c("HistSSB", "iSCAM_SSB0", "SSBMSY", "depletion", "50%Rmax", "90%RS")
frac <- c(1, 0.3, 0.4, 0.3, 1, 1)
HistSSB_y <- 2000

setup(length(type))
sfExportAll()
sfLibrary(dplyr)
sfLibrary(MSEtool)
FSH <- lapply(1:length(type), function(x) {
  Fopt(Hist, type = type[x], frac = frac[x], F_type = "abs", HistSSB_y = HistSSB_y, dep_type = "Dynamic")
})
saveRDS(FSH, file = "FSH_pcod.rds")


##### Second step - test two different HCRs

# Cutoff = SSB in 2000
# USR = mean SSB during 1956-2004
# Reference removal = mean F during 1956-2004
HCR_pcod <- function(Assessment, reps = 1, y_cutoff = 2000, y_USR = 1956:2004, y_RR = 1956:2004) {
  if(Assessment@conv) {
    SSB <- Assessment@SSB
    Year <- as.numeric(names(SSB))
    
    OCP_val <- SSB[length(SSB)]
    OCP <- c(mean(SSB[match(y_cutoff, Year)]), mean(SSB[match(y_USR, Year)]))
    FM <- Assessment@FMort
    RR <- mean(FM[match(y_RR, Year)])
    
    alpha <- SAMtool:::HCRlinesegment(OCP_val, OCP, relF = c(0, 1))
    TAC <- SAMtool:::calculate_TAC(Assessment, Ftarget = alpha * RR)
  } else {
    TAC <- NA_real_
  }
  Rec <- new("Rec")
  Rec@TAC <- rep(MSEtool::TACfilter(TAC), reps)
  return(Rec)
}
class(HCR_pcod) <- "HCR"

pcod_assess <- function(x, Data, MW = TRUE, start = list(R0 = 1.6e4, sigma_W = 0.15)) {
  if(max(Data@Year) > Data@LHYear) {
    Data@AddInd[, 1:4, match(Data@LHYear + 1, Data@Year):length(Data@Year)] <- NA_real_
    
    if(MW) {
      MW_sim <- apply(Data@CAL[x, , ][match(Data@LHYear + 1, Data@Year):length(Data@Year), , drop = FALSE], 1, function(xx) {
        weighted.mean(Data@wla[x]*Data@CAL_mids^Data@wlb[x], xx, na.rm = TRUE)
      })
      age <- 0:Data@MaxAge
      a <- Data@wla[x]
      b <- Data@wlb[x]
      Linf <- Data@vbLinf[x]
      K <- Data@vbK[x]
      t0 <- Data@vbt0[x]
      La <- Linf * (1 - exp(-K * (age - t0)))
      Wa <- a * La ^ b
      MW_sim <- apply(Data@CAA[x, , ][match(Data@LHYear + 1, Data@Year):length(Data@Year), , drop = FALSE], 1, function(xx) {
        weighted.mean(Wa, xx, na.rm = TRUE)
      })
      
      Data@Misc[[x]]$MW <- c(pcod$data@MS[, 1], MW_sim)
    }
  } else if(MW) {
    Data@Misc[[x]]$MW <- pcod$data@MS[, 1]
  }
  DD_SS(x, Data, MW = MW, start = start, prior = pcod$prior, AddInd = 1:5)
}
class(pcod_assess) <- "Assess"

#setup(8)
#prelim_AM(Hist@Data, pcod_assess)
#prelim_AM(Hist@Data, SCA)

pcod_MP_perfect <- make_MP(Perfect, HCR_pcod)
pcod_MP_DD <- make_MP(pcod_assess, HCR_pcod, MW = FALSE)
#pcod_MP_shortcut <- make_MP(Shortcut2, HCR_pcod)

Hist@OM@interval <- 2
#debug(HCR_pcod)
sfExport(list = c("pcod_assess", "HCR_pcod"))
MSE2 <- Project(Hist, MPs = c("pcod_MP_DD", "pcod_MP_perfect"), parallel = FALSE, extended = TRUE, checkMPs = FALSE)
saveRDS(MSE2, file = "Pcod/Pcod_HCR2.rds")

MSE2 <- readRDS("Pcod/Pcod_HCR.rds")
diagnostic(MSE2)

# Probability of being above 5 definitions of serious harm 
type <- c("HistSSB", "SSBMSY", "depletion", "50%Rmax", "90%RS")
frac <- c(1, 0.4, 0.3, 1, 1)
HistSSB_y <- 2000


PM <- Map(function(x, y) SHPM(MSE2, type = x, frac = y, HistSSB_y = HistSSB_y), x = type, y = frac)

annual_PM <- lapply(PM, function(x) apply(x@Stat >= x@Ref, 2:3, mean) %>% 
                      structure(dimnames = list(MP = MSE2@MPs, Year = 2020 + 1:MSE2@proyears)) %>% reshape2::melt() %>% 
                      dplyr::mutate(Metric = x@Caption)) %>%
  dplyr::bind_rows()

ggplot(annual_PM, aes(Year, value, colour = Metric)) + geom_point() + geom_line() + theme_bw() + 
  facet_wrap(~ MP) + 
  scale_colour_discrete(labels = scales::label_parse()) +
  labs(y = "Probability above serious harm", colour = "Serious harm candidates")
ggsave("Pcod/Pcod_SH_MP.png", height = 4, width = 10)

#terminal_PM <- dplyr::filter(annual_PM, Year == max(Year))

#ggplot(terminal_PM, aes(FM, value, colour = Metric)) + geom_point() + geom_line() + theme_bw() + 
#  scale_colour_discrete(labels = scales::label_parse()) +
#  labs(x = expression(F/F[MSY]), y = "Probability above serious harm in 2070", colour = "Serious harm candidates")
#ggsave("Pcod_SH_terminal.png", height = 4, width = 8)
