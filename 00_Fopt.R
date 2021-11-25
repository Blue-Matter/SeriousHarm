
library(MSEtool)
library(RPC)

optFSH <- function(logrelF, StockPars, FleetPars, Ncurr, nyears, proyears, nsim, F_type, p, SSBvec, opt = TRUE) {
  
  # optimize F, F/FMSY, F/M that provides probability p that stock > LRP
  if(F_type == "FMSY") {
    FM <- exp(logrelF) * StockPars$FMSY_y[, nyears + proyears]
  } else if(F_type == "M") {
    FM <- exp(logrelF) * StockPars$Marray[, nyears + proyears]
  } else {
    FM <- exp(logrelF) %>% rep(nsim)
  }
  
  simpop <- lapply(1:nsim, function(x) {
    MSEtool:::popdynCPP(nareas = StockPars$nareas, 
              maxage = StockPars$maxage, 
              Ncurr = Ncurr[x, , ], 
              pyears = proyears,
              M_age = StockPars$M_ageArray[x,,(nyears):(nyears+proyears)],
              Asize_c = StockPars$Asize[x,],
              MatAge=StockPars$Mat_age[x,,(nyears):(nyears+proyears)],
              WtAge=StockPars$Wt_age[x,,(nyears):(nyears+proyears)],
              FecAge=StockPars$Fec_Age[x,,(nyears):(nyears+proyears)],
              Vuln=FleetPars$V_real[x,,(nyears):(nyears+proyears)],
              Retc=FleetPars$retA_real[x,,(nyears):(nyears+proyears)],
              Prec=StockPars$Perr_y[x,(nyears):(nyears+proyears+StockPars$maxage)],
              movc=MSEtool:::split.along.dim(StockPars$mov[x,,,,(nyears):(nyears+proyears)],4),
              SRrelc=StockPars$SRrel[x],
              Effind=FleetPars$Find[x,],
              Spat_targc=FleetPars$Spat_targ[x],
              hc=StockPars$hs[x],
              R0c=StockPars$R0a[x,],
              SSBpRc=StockPars$SSBpR[x,],
              aRc=StockPars$aR[x,],
              bRc=StockPars$bR[x,],
              Qc=0,
              Fapic=FM[x],
              MPA=FleetPars$MPA,
              maxF=StockPars$maxF,
              control=2,
              SSB0c=StockPars$SSB0[x],
              plusgroup=StockPars$plusgroup)
  })
  SSB <- sapply(simpop, function(y) y[[4]], simplify = "array") %>% apply(c(4, 2), sum)
  N <- sapply(simpop, function(y) y[[1]], simplify = "array") %>% apply(c(4, 1, 2), sum)
  
  if(opt) {
    return((mean(SSB[, proyears] >= SSBvec) - p)^2)
  } else {
    return(list(SSB = SSB, N = N))
  }
}


Fopt <- function(Hist, type = c("HistSSB", "SSBMSY", "depletion", "50%Rmax", "90%RS", "iSCAM_SSB0"), frac = 1, 
                 F_type = c("FMSY", "FM", "abs"), p = 0.95, ...) {
  type <- match.arg(type)
  F_type <- match.arg(F_type)
  dots <- list(...)
  
  SSBvec <- do.call(get_SSBref, c(list(Hist = Hist, type = type), dots))
  StockPars <- Hist@SampPars$Stock
  StockPars$FMSY_y <- Hist@Ref$ByYear$FMSY
  FSH <- optimize(optFSH, interval = log(c(1e-8, 10)), StockPars = StockPars, FleetPars = Hist@SampPars$Fleet, 
                  Ncurr = Hist@AtAge$Number[, , Hist@OM@nyears, ], 
                  nyears = Hist@OM@nyears, proyears = Hist@OM@proyears, nsim = Hist@OM@nsim,
                  F_type = F_type, p = p, SSBvec = frac * SSBvec)
  
  out <- optFSH(logrelF = exp(FSH$minimum), StockPars = StockPars, FleetPars = Hist@SampPars$Fleet, 
                Ncurr = Hist@AtAge$Number[, , Hist@OM@nyears, ], 
                nyears = Hist@OM@nyears, proyears = Hist@OM@proyears, nsim = Hist@OM@nsim,
                F_type = F_type, p = p, SSBvec = frac * SSBvec, opt = FALSE)

  c(list(opt = FSH), out)
}


get_SSBref <- function(Hist, type = c("HistSSB", "SSBMSY", "depletion", "50%Rmax", "90%RS", "iSCAM_SSB0"), ...) {
  type <- match.arg(type)
  dots <- list(...)
  
  OM <- Hist@OM
  
  yr_pro <- seq(OM@CurrentYr + 1, OM@CurrentYr + OM@proyears)
  #yind <- match(dots$year_range[1]:dots$year_range[2], yr)
  
  
  if(type == "HistSSB") {
    yr_hist <- seq(OM@CurrentYr - OM@nyears + 1, OM@CurrentYr)
    xout <- Hist@TSdata$SBiomass %>% apply(1:2, sum)
    ref <- xout[, match(dots$HistSSB_y, yr_hist)]
  } else if(type == "SSBMSY") {
    xout <- Hist@Ref$ByYear$SSBMSY
    ref <- xout[, ncol(xout)]
  } else if(type == "depletion") {
    
    ref <- switch(dots$dep_type,
                  "Equilibrium" = Hist@Ref$ByYear$SSB0[, OM@nyears + OM@proyears],
                  "Initial" = Hist@Ref$ByYear$SSB0[, 1],
                  "Dynamic" = Hist@Ref$Dynamic_Unfished$SSB0[, OM@nyears + OM@proyears])
  } else if(type == "50%Rmax") {
    ref <- RPC:::calculate_SSB50(Hist)$SSB50
  } else if(type == "90%RS") {
    ref <- local({
      out <- RPC:::stock_recruit_int(Hist)
      RpS_90 <- apply(out$R/out$SSB, 1, quantile, probs = 0.9)
      R_90 <- apply(out$R, 1, quantile, probs = 0.9)
      R_90/RpS_90
    })
  } else {
    
    .iSCAM_SSB0 <- function(Hist) {
      nyears <- Hist@OM@nyears
      maxage <- Hist@OM@maxage
      
      M <- apply(Hist@SampPars$Stock$M_ageArray[, , 1:nyears], 1:2, mean)
      Wt <- apply(Hist@SampPars$Stock$Wt_age[, , 1:nyears], 1:2, mean)
      Mat <- apply(Hist@SampPars$Stock$Mat_age[, , 1:nyears], 1:2, mean)
      Fec <- apply(Hist@SampPars$Stock$Fec_Age[, , 1:nyears], 1:2, mean)
      V <- rep(1, maxage + 1)
      R0 <- Hist@SampPars$Stock$R0
      h <- Hist@SampPars$Stock$hs
      SRrel <- Hist@OM@SRrel
      phi <- Hist@SampPars$Stock$SSBpR[, 1]
      
      sapply(1:Hist@OM@nsim, function(x) {
        MSEtool:::MSYCalcs(logF = 0, M_at_Age = M[x, ], Wt_at_Age = Wt[x, ], Mat_at_Age = Mat[x, ], Fec_at_Age = Fec[x, ], V_at_Age = V, 
                           maxage = maxage, R0x = R0[x], SRrelx = SRrel, hx = h[x], SSBpR = phi[x], opt = 0L)["SB0"]
      })
    }
    ref <- .iSCAM_SSB0(Hist)
  }
  return(ref)
}

