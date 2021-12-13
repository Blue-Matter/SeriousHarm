

library(SAMtool)

nsim <- 48
OM <- RPC::DFO_WCVI_Herring_2019 %>% SubCpars(1:nsim)
OM@cpars$Data@vbt0 <- NA
OM@interval <- 200
Hist <- runMSE(OM, Hist = TRUE, parallel = TRUE)


##### First step - fixed F MPs to identify suitability of serious harm metrics to fishing
MPs <- paste0(relF, "_FMSY")
MSE <- Project(Hist, MPs = MPs, extended = TRUE)
saveRDS(MSE, file = "Pherring/Pherring_FMSY.rds")

MSE <- readRDS("Pherring/Pherring_FMSY.rds")

##### First step - fixed F MPs to identify suitability of serious harm metrics to fishing
MPs <- paste0(relF, "_M")
MSE <- Project(Hist, MPs = MPs, extended = TRUE)
saveRDS(MSE, file = "Pherring/Pherring_M.rds")

MSE <- readRDS("Pherring/Pherring_M.rds")


##### Calculate F such that P(SH) = 0.95
source("00_Fopt.R")

# Show historical probability of serious harm
type <- c("HistSSB", "iSCAM_SSB0", "SSBMSY", "depletion", "50%Rmax", "90%RS")
frac <- c(1, 0.3, 0.4, 0.3, 1, 1)
HistSSB_y <- 2010

if(!snowfall::sfIsRunning()) setup(length(type))
sfExportAll()
sfLibrary(dplyr)
sfLibrary(MSEtool)
FSH_abs <- sfLapply(1:length(type), function(x) {
  Fopt(Hist, LRP_type = type[x], frac = frac[x], F_type = "abs", HistSSB_y = HistSSB_y, dep_type = "Dynamic", p = 0.95)
})
saveRDS(FSH_abs, file = "Pherring/Pherring_FP95_abs.rds")
FSH_M <- sfLapply(1:length(type), function(x) {
  Fopt(Hist, LRP_type = type[x], frac = frac[x], F_type = "M", HistSSB_y = HistSSB_y, dep_type = "Dynamic", p = 0.95, F_search = seq(0.01, 4, 0.01))
})
saveRDS(FSH_M, file = "Pherring/Pherring_FP95_M.rds")
sfStop()

##### Second step - test two different HCRs
# See Table 1 of https://waves-vagues.dfo-mpo.gc.ca/Library/40760133.pdf
HCR_pherring <- function(Assessment, reps = 1, HCR_type = c("escapement", "hockeystick"), OCP_type = c("SSB", "depletion"),
                         SSB_cutoff = 18.8, SSB_upper = NULL, max_harvest = 0.2, max_cap = NULL) {
  HCR_type <- match.arg(HCR_type)
  OCP_type <- match.arg(OCP_type)
  
  if(Assessment@conv) {
    SSB <- Assessment@SSB
    Year <- as.numeric(names(SSB))
    
    if(OCP_type == "SSB") {
      OCP_val <- SSB[length(SSB)]
    } else {
      SSB0 <- local({
        nhist <- 69
        M <- apply(Assessment@TMB_report$M[1:nhist, ], 2, mean)
        wt <- Assessment@info$data$weight
        mat <- Assessment@info$data$mat
        NPR <- numeric(length(wt))
        NPR[1] <- 1
        for(a in 2:length(NPR)) NPR[a] <- NPR[a-1] * exp(-M[a-1])
        NPR[length(NPR)] <- NPR[length(NPR)]/(1 - exp(-M[length(NPR)]))
        phi0 <- sum(NPR * mat * wt)
        
        alpha <- Assessment@TMB_report$Arec
        beta <- Assessment@TMB_report$Brec
        
        num <- alpha * phi0 - 1 # BH only
        num/beta
      })
      
      OCP_val <- SSB[length(SSB)]/SSB0
    }
    
    if(HCR_type == "escapement") {
      OCP <- rep(SSB_cutoff, 2)
    } else {
      OCP <- c(SSB_cutoff, SSB_upper)
    }
    F_harvest <- -log(1 - max_harvest)
    relF <- c(0, F_harvest)
    
    alpha <- SAMtool:::HCRlinesegment(OCP_val, OCP, relF)
    TAC <- SAMtool:::calculate_TAC(Assessment, Ftarget = alpha * F_harvest)
    if(is.null(max_cap)) TAC <- min(TAC, max_cap)
  } else {
    TAC <- NA_real_
  }
  Rec <- new("Rec")
  Rec@TAC <- rep(MSEtool::TACfilter(TAC), reps)
  return(Rec)
}
class(HCR_pherring) <- "HCR"

pherring_assess <- function(x, Data, ...) {
  if(max(Data@Year) > Data@LHYear) {
    Data@AddInd[, 1, match(Data@LHYear, Data@Year):length(Data@Year)] <- NA_real_
  }
  SAMtool::SCA_RWM(x, Data, prior = list(M = c(Data@Misc$StockPars$M_ageArray[x, 1, 1], 0.2)), AddInd = 1:2, M_bounds = c(0.1, 2))
}
class(pherring_assess) <- "Assess"

MP1 <- make_MP(pherring_assess, HCR_pherring, HCR_type = "escapement", OCP_type = "SSB", max_harvest = 0.2)
MP8 <- make_MP(pherring_assess, HCR_pherring, HCR_type = "hockeystick", OCP_type = "depletion",
               SSB_cutoff = 0.3, SSB_upper = 0.6, max_harvest = 0.1, max_cap = NULL)

setup(12)
sfExport(list = c("pherring_assess", "HCR_pherring", "MP1", "MP8"))
Hist@OM@interval <- 2
MSE2 <- Project(Hist, MPs = c("MP1", "MP8"), extended = TRUE, parallel = TRUE, checkMPs = FALSE)
saveRDS(MSE2, file = "Pherring/Pherring_HCR.rds")

MSE2 <- readRDS("Pherring/Pherring_HCR.rds")
diagnostic(MSE2)

# Probability of being above 5 definitions of serious harm
type <- c("HistSSB", "iSCAM_SSB0", "SSBMSY", "depletion", "50%Rmax", "90%RS")
frac <- c(1, 0.3, 0.4, 0.4, 1, 1)

PM <- Map(function(x, y) SHPM(MSE2, type = x, frac = y, HistSSB_y = 2010, dep_type = "Dynamic"), x = type, y = frac)

annual_PM <- lapply(PM, function(x) apply(x@Stat >= x@Ref, 2:3, mean) %>% 
                      structure(dimnames = list(MP = MSE2@MPs, Year = MSE2@OM$CurrentYr[1] + 1:MSE2@proyears)) %>% reshape2::melt() %>% 
                      dplyr::mutate(Metric = x@Caption)) %>%
  dplyr::bind_rows()

ggplot(annual_PM, aes(Year, value, colour = Metric)) + geom_point() + geom_line() + theme_bw() + 
  facet_wrap(~ MP) + 
  scale_colour_discrete(labels = scales::label_parse()) +
  labs(y = "Probability above serious harm", colour = "Serious harm candidates") + coord_cartesian(ylim = c(0, 1))
ggsave("Pherring/Pherring_SH_MP.png", height = 4, width = 10)

#terminal_PM <- dplyr::filter(annual_PM, Year == max(Year))

#ggplot(terminal_PM, aes(FM, value, colour = Metric)) + geom_point() + geom_line() + theme_bw() + 
#  scale_colour_discrete(labels = scales::label_parse()) +
#  labs(x = expression(F/F[MSY]), y = "Probability above serious harm in 2070", colour = "Serious harm candidates")
#ggsave("Pcod_SH_terminal.png", height = 4, width = 8)
