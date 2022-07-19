

library(MSEtool)
source("analysis_fn.R")


# SWO 
OM <- readRDS(file = "OM/OM_ICCAT_SWO.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

SW <- get_ref_pt(Hist, OM, 
                 Year_assess = 1999,
                 Name = "N.A.~Swordfish",
                 Yr_SRP = 2002,  # Year corresponding to ICES Blim (Type 5, no B_SR)
                 Yr_Fmed = 1979, # First year for calculating median R/S
)
#SW$RPts$SB_recover <- SW$Assess$SB_recover <- NA
saveRDS(SW, file = "LRP/LRP_ICCAT_SWO.rds")
SW <- readRDS(file = "LRP/LRP_ICCAT_SWO.rds")

plot_F(SW)
plot_B(SW)
plot_SR(SW)
plot_Hist(SW)
plot_SR_LRP(SW)
plot_SP(SW)

# RS
OM <- readRDS(file = "OM/OM_GoM_RS.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

RS <- get_ref_pt(Hist, OM, 
                 Year_assess = 1990, # 1990 Gulf FMP Amendment 1, 1991 Amendment 3
                 Name = "GoM~Red~Snapper",
                 Yr_SRP = "hockey_stick",  # Year corresponding to ICES Blim (Type 4, no B_SR)
                 Yr_Fmed = 1972 # First year for calculating median R/S
)

saveRDS(RS, file = "LRP/LRP_GoM_RS.rds")
RS <- readRDS(file = "LRP/LRP_GoM_RS.rds")

plot_F(RS)
plot_B(RS)
plot_SR(RS)
plot_Hist(RS)
plot_SR_LRP(RS)
plot_SP(RS)

# NE Acadian redfish
OM <- readRDS(file = "OM/OM_NE_redfish.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

AR <- get_ref_pt(Hist, OM, 
                 Year_assess = 1976, 
                 Name = "N.E.~Acadian~redfish",
                 Yr_SRP = 1971,  # Year corresponding to ICES Blim
                 Yr_Fmed = 1960,  # First year for calculating median R/S
                 Yr_SP = 1974
)


saveRDS(AR, file = "LRP/LRP_NE_redfish.rds")
AR <- readRDS(file = "LRP/LRP_NE_redfish.rds")

plot_F(AR)
plot_B(AR)
plot_SR(AR)
plot_Hist(AR)
plot_SR_LRP(AR)
plot_SP(AR)

# GoM Haddock
#OM <- readRDS("OM/OM_GoM_haddock.rds")
#Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)
#
## No formal assessment between 1986-2001 overfished statement 
## See Amendment 9 (1999)
#HD <- get_ref_pt(Hist, OM, 
#                 Year_assess = 1999, 
#                 Name = "G.M.~haddock",
#                 Yr_SRP = 2010,  # Year corresponding to ICES Blim
#                 #Yr_Fmed = 1980  # First year for calculating median R/S
#)
#saveRDS(HD, file = "LRP/LRP_GoM_haddock.rds")
#
#plot_F(HD)
#plot_B(HD)
#plot_SR(HD)
#plot_Hist(HD)
#plot_SR_LRP(HD)
#plot_SP(HD)

# Southern New England yellowtail flounder
OM <- readRDS("OM/OM_sne_yt.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)


SNEYT <- get_ref_pt(Hist, OM, 
                    Year_assess = 2005, # Framework 42
                    Name = "S.N.E.~yellowtail~flounder",
                    Yr_SRP = 1987,  # Year corresponding to ICES Blim
                    Yr_SP = 1987,
                    #Yr_Fmed = 1980  # First year for calculating median R/S
)
# B recover is SSB in 1987
SNEYT$RPvars$Yr_Brecover <- 1987
SNEYT$RPts$SB_recover <- SNEYT$RPts$SB[SNEYT$RPts$Year == 1987]
SNEYT$Assess$SB_recover <- SNEYT$RPts$SB[SNEYT$RPts$Year == 2005]/SNEYT$RPts$SB[SNEYT$RPts$Year == 1987]

saveRDS(SNEYT, file = "LRP/LRP_sne_yt.rds")

plot_F(SNEYT)
plot_B(SNEYT)
plot_SR(SNEYT)
plot_Hist(SNEYT)
plot_SR_LRP(SNEYT)
plot_SP(SNEYT)



# SA red porgy
OM <- readRDS(file = "OM/OM_SA_RedPorgy.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

RP <- get_ref_pt(Hist, OM, 
                 Year_assess = 1992, # Amendment 4 (see SEDAR 60 report)
                 Name = "S.A.~Red~Porgy",
                 Yr_SRP = "hockey_stick",  # Year corresponding to ICES Blim
                 #Yr_Fmed = 1990 # First year for calculating median R/S
)
RP$Assess$SB_recover <- RP$RPts$SB_recover <- RP$RPvars$Yr_Brecover <- NA

saveRDS(RP, file = "LRP/LRP_SA_RedPorgy.rds")
RP <- readRDS(file = "LRP/LRP_SA_RedPorgy.rds")

plot_F(RP)
plot_B(RP)
plot_SR(RP)
plot_Hist(RP)
plot_SR_LRP(RP, year = TRUE)
plot_SP(RP)


# Cowcod
OM <- readRDS(file = "OM/OM_SCB_cowcod.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

CC <- get_ref_pt(Hist, OM, 
                  Year_assess = 2000, # Amendment 4 (see SEDAR 60 report)
                  Name = "S.C.B.~Cowcod",
                  Yr_SRP = NA #,  # Year corresponding to ICES Blim
                  #Yr_Fmed = 1983 # First year for calculating median R/S
)
CC$Assess$SB_SRP <- CC$RPts$SB_SRP <- NA
saveRDS(CC, file = "LRP/LRP_SCB_cowcod.rds")
CC <- readRDS(file = "LRP/LRP_SCB_cowcod.rds")

plot_F(CC)
plot_B(CC)
plot_SR(CC)
plot_Hist(CC)

plot_SR_LRP(CC)
plot_SP(CC)


# NEA horse mackerel
OM <- readRDS(file = "OM/OM_NEA_hom.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

HOM <- get_ref_pt(Hist, OM, 
                  Year_assess = 2020, #?
                  Name = "N.E.A.~Horse~Mackerel",
                  Yr_SRP = 2017,  # Year corresponding to ICES Blim
                  Yr_Fmed = 1982 # First year for calculating median R/S
)
HOM$Assess$SB_recover <- HOM$RPts$SB_recover <- HOM$RPvars$Yr_Brecover <- NA
saveRDS(HOM, file = "LRP/LRP_NEA_hom.rds")
HOM <- readRDS(file = "LRP/LRP_NEA_hom.rds")

plot_F(HOM)
plot_B(HOM)

local({
  HOM$RPts <- filter(HOM$RPts, Year >= 1983)
  plot_SR_LRP(HOM)
})

plot_SR(HOM)
plot_SR_LRP(HOM)
plot_Hist(HOM)

plot_SP(HOM)


# US Dark-blotched rockfish
OM <- readRDS(file = "OM/OM_WC_darkblotched.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)


DB <- get_ref_pt(Hist, OM, 
                 Year_assess = 2003, # Declared overfished
                 Name = "U.S.~Darkblotched~rockfish",
                 Yr_SRP = 2010,  # Year corresponding to ICES Blim
                 Yr_Fmed = 1970 # First year for calculating median R/S
)
DB$Assess$SB_recover <- DB$RPts$SB_recover <- DB$RPvars$Yr_Brecover <- NA

saveRDS(DB, file = "LRP/LRP_WC_darkblotched.rds")
DB <- readRDS(file = "LRP/LRP_WC_darkblotched.rds")


plot_F(DB)
plot_B(DB)
plot_SR(DB)
plot_Hist(DB)
plot_SR_LRP(DB)

plot_SP(DB)

# US POP
OM <- readRDS(file = "OM/OM_WC_pop.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

POP <- get_ref_pt(Hist, OM, 
                  Year_assess = 2001, # Declared overfished
                  Name = "U.S.~Pacific~ocean~perch",
                  Yr_SRP = 1987,  # Year corresponding to ICES Blim
                  Yr_Fmed = 1965 # First year for calculating median R/S
)

saveRDS(POP, file = "LRP/LRP_WC_pop.rds")
POP <- readRDS(file = "LRP/LRP_WC_pop.rds")

local({
  POP$RPts <- filter(POP$RPts, Year != 2008)
  plot_SR_LRP(POP)
})

plot_F(POP)
plot_B(POP)
plot_SR(POP)
plot_Hist(POP)

plot_SR_LRP(POP)
plot_SP(POP)


# sGSL cod
OM <- readRDS(file = "OM/OM_sGSL_cod.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

GSLC <- get_ref_pt(Hist, OM, 
                   Year_assess = 1993, # See 2005 SAR
                   Name = "s.GSL~Atl.~Cod",
                   Yr_SRP = "hockey_stick",  # Year corresponding to ICES Blim
                   Yr_Fmed = 1971, # First year for calculating median R/S,
                   Yr_SP = 1992
)

#GSLC$RPts$FMSY <- GSLC$RPts$SB0 <- GSLC$RPts$
saveRDS(GSLC, file = "LRP/LRP_sGSL_cod.rds")
GSLC <- readRDS(file = "LRP/LRP_sGSL_cod.rds")


plot_F(GSLC)
plot_B(GSLC)
plot_SR(GSLC)
plot_Hist(GSLC)
plot_SR_LRP(GSLC, year = TRUE)
plot_SP(GSLC)

# sGSL herring
OM <- readRDS(file = "OM/OM_sGSL_herring.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

# Update hist for change in SR rel
rp <- sapply(1:OM@nyears, function(y) {
  
  opt2 <- optimize(SAMtool:::yield_fn_SCA, interval = c(1e-4, 4), 
                   M = Hist@SampPars$Stock$M_ageArray[1, , y], 
                   mat = Hist@SampPars$Stock$Mat_age[1, , y], 
                   weight = Hist@SampPars$Stock$Wt_age[1, , y], 
                   vul = Hist@SampPars$Fleet$V_real[1, , y], 
                   SR = "Ricker", 
                   Arec = ifelse(y >= 15, 8.48, 19.89), 
                   Brec = ifelse(y >= 15, 9.7e-6, 1.21e-5), 
                   catch_eq = "Baranov")
  
  opt3 <- SAMtool:::yield_fn_SCA(opt2$minimum,
                                 M = Hist@SampPars$Stock$M_ageArray[1, , y], 
                                 mat = Hist@SampPars$Stock$Mat_age[1, , y], 
                                 weight = Hist@SampPars$Stock$Wt_age[1, , y], 
                                 vul = Hist@SampPars$Fleet$V_real[1, , y], 
                                 SR = "Ricker", 
                                 Arec = ifelse(y >= 15, 8.48, 19.89), 
                                 Brec = ifelse(y >= 15, 9.7e-6, 1.21e-5), 
                                 catch_eq = "Baranov",
                                 opt = FALSE)
  FMSY <- opt2$minimum
  MSY <- -1 * opt2$objective
  SBMSY <- opt3["E"]
  
  unfished <- SAMtool:::yield_fn_SCA(0,
                                     M = Hist@SampPars$Stock$M_ageArray[1, , y], 
                                     mat = Hist@SampPars$Stock$Mat_age[1, , y], 
                                     weight = Hist@SampPars$Stock$Wt_age[1, , y], 
                                     vul = Hist@SampPars$Fleet$V_real[1, , y], 
                                     SR = "Ricker", 
                                     Arec = ifelse(y >= 15, 8.48, 19.89), 
                                     Brec = ifelse(y >= 15, 9.7e-6, 1.21e-5), 
                                     catch_eq = "Baranov",
                                     opt = FALSE)
  SB0 <- unfished["E"]
  
  c("FMSY" = FMSY, "SBMSY" = SBMSY, "SB0" = SB0)
})

Hist@Ref$ByYear$FMSY[1, 1:OM@nyears] <- rp[1, ]
Hist@Ref$ByYear$SSBMSY[1, 1:OM@nyears] <- rp[2, ]
Hist@Ref$ByYear$SSB0[1, 1:OM@nyears] <- rp[3, ]

GSLH <- get_ref_pt(Hist, OM, 
                   Year_assess = 2018, #?
                   Name = "s.GSL~Atl.~Herring",
                   Yr_SRP = NA,  # Year corresponding to ICES Blim
                   #Yr_Fmed = 1983 # First year for calculating median R/S
)


saveRDS(GSLH, file = "LRP/LRP_sGSL_herring.rds")
GSLH <- readRDS(file = "LRP/LRP_sGSL_herring.rds")


plot_F(GSLH)
plot_B(GSLH)
plot_SR(GSLH)
plot_Hist(GSLH)

plot_SR_LRP(GSLH)
plot_SP(GSLH)


# NAFO plaice
OM <- readRDS(file = "OM/OM_NAFO_plaice.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

# Calculate reference points with M = 0.2 in years where population M = 0.53 (temporary)
for(y in 1:OM@nyears) {
  if(any(OM@cpars$M_ageArray[1, , y] > 0.2)) {
    
    MSYs <- local({
      OM@cpars$M_ageArray[1, , y] <- OM@cpars$M_ageArray[1, , 1]
      optMSY_eq(1, 
                M_ageArray = OM@cpars$M_ageArray,
                Wt_age = OM@cpars$Wt_age,
                Mat_age = OM@cpars$Mat_age,
                Fec_age = OM@cpars$Mat_age * OM@cpars$Wt_age,
                V = OM@cpars$V,
                maxage = OM@maxage,
                R0 = OM@cpars$R0,
                SRrel = OM@SRrel,
                hs = OM@cpars$hs,
                SSBpR = Hist@SampPars$Stock$SSBpR,
                yr.ind = y,
                plusgroup = 1)
    })
    
    Fcalc <- local({
      OM@cpars$M_ageArray[1, , y] <- OM@cpars$M_ageArray[1, , 1]
      boundsF <- c(1E-3, 3)
      F_search <- exp(seq(log(min(boundsF)), log(max(boundsF)), length.out = 50))
      Ref_search <- MSEtool:::Ref_int_cpp(F_search, 
                                M_at_Age = OM@cpars$M_ageArray[1, , y],
                                Wt_at_Age = OM@cpars$Wt_age[1, , y], 
                                Mat_at_Age = OM@cpars$Mat_age[1, , y],
                                Fec_at_Age = OM@cpars$Wt_age[1, , y] * OM@cpars$Mat_age[1, , y],
                                V_at_Age = OM@cpars$V[1, , y],
                                maxage = OM@maxage,
                                plusgroup = 1)
      RPS <- Ref_search[3,]
      SSB <- Hist@TSdata$SBiomass[1, , ] %>% rowSums()
      R <- Hist@AtAge$Number[1, 1, , ] %>% rowSums()
      Fmed <- MSEtool:::LinInterp_cpp(RPS, F_search, xlev = median(R/SSB))
      Fmed
    })
    
    Hist@Ref$ByYear$FMSY[1, y] <- MSYs["F"]
    Hist@Ref$ByYear$SSBMSY[1, y] <- MSYs["SB"]
    Hist@Ref$ByYear$SSB0[1, y] <- MSYs["SB0"]
    
    Hist@Ref$ByYear$Fmed[1, y] <- Fcalc
  }
}

PL <- get_ref_pt(Hist, OM, 
                 Year_assess = 1994, # Fishing moratorium
                 Name = "NAFO~American~plaice",
                 Yr_SRP = "hockey_stick",  # Year corresponding to ICES Blim
                 #Yr_Fmed = 1960, # First year for calculating median R/S
                 Yr_SP = 1993
)


PL$RPts$SB_recover <- PL$RPts$SB[17]
PL$Assess$SB_recover <- PL$RPts$SB[35]/ PL$RPts$SB[17]
PL$RPvars$Yr_Brecover <- 1976

saveRDS(PL, file = "LRP/LRP_NAFO_plaice.rds")
PL <- readRDS(file = "LRP/LRP_NAFO_plaice.rds")


plot_F(PL)
plot_B(PL)
plot_SR(PL)
ggsave("SR_plaice.png", height = 4, width = 5)

plot_Hist(PL)
plot_SR_LRP(PL)

plot_SP(PL)



# NAFO cod
OM <- readRDS(file = "OM/OM_NAFO_cod.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

NAFOcod <- get_ref_pt(Hist, OM, 
                      Year_assess = 1994, # Fishing moratorium
                      Name = "NAFO~Atlantic~cod",
                      Yr_SRP = "hockey_stick",  # Year corresponding to ICES Blim
                      #Yr_Fmed = 1960, # First year for calculating median R/S
                      Yr_SP = 1991
)
NAFOcod$RPts$SB_recover <- NAFOcod$RPts$SB[18]
NAFOcod$Assess$SB_recover <- NAFOcod$RPts$SB[36]/ NAFOcod$RPts$SB[18]
NAFOcod$RPvars$Yr_Brecover <- 1976

saveRDS(NAFOcod, file = "LRP/LRP_NAFO_cod.rds")
NAFOcod <- readRDS(file = "LRP/LRP_NAFO_cod.rds")


plot_F(NAFOcod)
plot_B(NAFOcod)
plot_SR(NAFOcod)
#ggsave("SR_plaice.png", height = 4, width = 5)

plot_Hist(NAFOcod)
plot_SR_LRP(NAFOcod)

plot_SP(NAFOcod)



# GoM cod
OM <- readRDS(file = "OM/OM_GoM_cod.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

GoM_cod <- get_ref_pt(Hist, OM, 
                      Year_assess = 2005, # Framework 42
                      Name = "G.M.~Atlantic~cod",
                      Yr_SRP = "hockey_stick",  # Year corresponding to ICES Blim
                      #Yr_Fmed = 1960, # First year for calculating median R/S
                      Yr_SP = 2011
)
GoM_cod$RPts$SB_recover <- GoM_cod$RPts$SB[13]
GoM_cod$Assess$SB_recover <- GoM_cod$RPts$SB[24]/ GoM_cod$RPts$SB[13]
GoM_cod$RPvars$Yr_Brecover <- 1994

saveRDS(GoM_cod, file = "LRP/LRP_GoM_cod.rds")
GoM_cod <- readRDS(file = "LRP/LRP_GoM_cod.rds")


plot_F(GoM_cod)
plot_B(GoM_cod)
plot_SR(GoM_cod)
plot_Hist(GoM_cod)
plot_SR_LRP(GoM_cod)

plot_SP(GoM_cod)


# BET
OM <- readRDS(file = "OM/OM_ICCAT_BET.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

BET <- get_ref_pt(Hist, OM, 
                  Year_assess = 2018, # ICCAT resolutions and lobbying
                  Name = "Atlantic~bigeye~tuna",
                  Yr_SRP = 2011,  # Year corresponding to ICES Blim
                  Yr_Fmed = 1973 # First year for calculating median R/S
)
BET$Assess$SB_SRP <- BET$RPts$SB_SRP <- NA

saveRDS(BET, file = "LRP/LRP_ICCAT_BET.rds")
BET <- readRDS(file = "LRP/LRP_ICCAT_BET.rds")


plot_F(BET)
plot_B(BET)
plot_SR(BET)
plot_Hist(BET)
plot_SR_LRP(BET)

plot_SP(BET)


# WCVI herring
OM <- readRDS(file = "OM/OM_WCVI_herring.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

WCVI <- get_ref_pt(Hist, OM, 
                   Year_assess = 2005, # close fishery
                   Name = "WCVI~Pacific~herring",
                   Yr_SRP = "hockey_stick",  # Year corresponding to ICES Blim
                   Yr_Fmed = 1951, # First year for calculating median R/S
                   Yr_SP = 2005
)

saveRDS(WCVI, file = "LRP/LRP_WCVI_herring.rds")
plot_F(WCVI)
plot_B(WCVI)
plot_SR(WCVI)
plot_Hist(WCVI)
plot_SR_LRP(WCVI)

plot_SP(WCVI)


# EBS cod
OM <- readRDS(file = "OM/OM_EBS_cod.rds")
Hist <- SubCpars(OM, 1:2) %>% runMSE(Hist = TRUE)

EBScod <- get_ref_pt(Hist, OM, 
                     Year_assess = 2019, # zero catch advised https://doi.org/10.17895/ices.advice.4747
                     Name = "E.B.S.~Atlantic~cod",
                     Yr_SRP = 2012,  # Year corresponding to ICES Blim
                     Yr_Fmed = 1946, # First year for calculating median R/S
                     Yr_SP = 2012
)
EBScod$RPts$SB_recover <- EBScod$RPts$SB[15]
EBScod$Assess$SB_recover <- EBScod$RPts$SB[74]/ EBScod$RPts$SB[15]
EBScod$RPvars$Yr_Brecover <- 1960


saveRDS(EBScod, file = "LRP/LRP_EBS_cod.rds")
plot_F(EBScod)
plot_B(EBScod)
plot_SR(EBScod)
plot_Hist(EBScod)
plot_SR_LRP(EBScod)

plot_SP(EBScod)

