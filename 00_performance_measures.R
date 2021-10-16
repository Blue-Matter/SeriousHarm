
library(MSEtool)

# Make functions for performance measures based on various definitions of serious harm
SHPM <- function(MSE, type = c("HistSSB", "SSBMSY", "depletion", "50%Rmax", "90%RS", "iSCAM_SSB0"), frac = 1, Hist = FALSE, ...) {
  type <- match.arg(type)
  dots <- list(...)
  
  if(is.null(dots$year_range)) {
    if(Hist) {
      dots$year_range <- c(MSE@OM$CurrentYr[1] - MSE@nyears + 1, MSE@OM$CurrentYr[1])
    } else {
      dots$year_range <- MSE@OM$CurrentYr[1] + c(1, MSE@proyears)
    }
  }
  #if(Hist) MSEhist <- MSE@Hist
  
  if(type == "HistSSB") {
    if(is.null(dots$HistSSB_y)) dots$HistSSB_y <- MSE@OM$CurrentYr[1]
    if(Hist) {
      PM <- local({
        yr <- seq(MSE@OM$CurrentYr[1] - MSE@nyears + 1, MSE@OM$CurrentYr[1])
        yind <- match(dots$year_range[1]:dots$year_range[2], yr)
        xout <- MSE@SSB_hist
        ref <- MSE@SSB_hist[, match(dots$HistSSB_y, yr)]
        caption <- paste0("Probability~SSB/SSB[", dots$HistSSB_y, "]>", frac)
        label <- caption %>% strsplit("~") %>% getElement(1) %>% paste(collapse = " ")
        
        PM <- new("PMobj")
        PM@Name <- label
        PM@Caption <- caption
        PM@Stat <- xout/ref
        PM@Ref <- frac
        PM@Prob <- apply(PM@Stat[, yind, drop = FALSE] > PM@Ref, 2, mean)
        PM@Mean <- mean(PM@Stat[, yind, drop = FALSE])
        PM@MPs <- "Hist"
        PM
      })
    } else {
      PM <- RPC::make_PMobj(x = MSE, type = "SSB", frac = frac, year_range = dots$year_range, SSBhist_yr = dots$HistSSB_y)
    }
    
  } else if(type == "SSBMSY") {
    
    if(Hist) {
      PM <- local({
        yr <- seq(MSE@OM$CurrentYr[1] - MSE@nyears + 1, MSE@OM$CurrentYr[1])
        yind <- match(dots$year_range[1]:dots$year_range[2], yr)
        xout <- MSE@SSB_hist
        ref <- MSE@RefPoint$ByYear$SSBMSY[, 1:MSE@nyears]
        caption <- paste0("Probability~SSB/SSB[MSY]>", frac)
        label <- caption %>% strsplit("~") %>% getElement(1) %>% paste(collapse = " ")
        
        PM <- new("PMobj")
        PM@Name <- label
        PM@Caption <- caption
        PM@Stat <- xout/ref
        PM@Ref <- frac
        PM@Prob <- apply(PM@Stat[, yind, drop = FALSE] > PM@Ref, 2, mean)
        PM@Mean <- mean(PM@Stat[, yind, drop = FALSE])
        PM@MPs <- "Hist"
        PM
      })
    } else {
      PM <- RPC::make_PMobj(x = MSE, type = "SSBMSY", frac = frac, year_range = dots$year_range)
    }
    
  } else if(type == "depletion") {
    
    if(is.null(dots$dep_type)) dots$dep_type <- "Equilibrium"
    
    if(Hist) {
      PM <- local({
        yr <- seq(MSE@OM$CurrentYr[1] - MSE@nyears + 1, MSE@OM$CurrentYr[1])
        yind <- match(dots$year_range[1]:dots$year_range[2], yr)
        xout <- MSE@SSB_hist
        ref <- switch(dots$dep_type,
                      "Equilibrium" = MSE@RefPoint$ByYear$SSB0[, 1:MSE@nyears],
                      "Dynamic" = MSE@RefPoint$Dynamic_Unfished$SSB0[, 1:MSE@nyears],
                      "Initial" = MSE@SSB_hist[, 1])
        caption <- paste0("Probability~'SSB'~'/'~'", dots$dep_type, "'~SSB[0]>", frac)
        
        label <- caption %>% strsplit("~") %>% getElement(1) %>% paste(collapse = " ")
        
        PM <- new("PMobj")
        PM@Name <- label
        PM@Caption <- caption
        PM@Stat <- xout/ref
        PM@Ref <- frac
        PM@Prob <- apply(PM@Stat[, yind, drop = FALSE] > PM@Ref, 2, mean)
        PM@Mean <- mean(PM@Stat[, yind, drop = FALSE])
        PM@MPs <- "Hist"
        PM
      })
    } else {
      PM <- RPC::make_PMobj(x = MSE, type = "SSB0", frac = frac, year_range = dots$year_range, SSB0_type = dots$dep_type)
    }
    
  } else if(type == "50%Rmax") {
    
    if(Hist) {
      PM <- local({
        yr <- seq(MSE@OM$CurrentYr[1] - MSE@nyears + 1, MSE@OM$CurrentYr[1])
        yind <- match(dots$year_range[1]:dots$year_range[2], yr)
        xout <- MSE@SSB_hist
        ref <- RPC:::calculate_SSB50(MSE@Hist)$SSB50
        caption <- paste0("Probability~SSB/SSB[\"50%Rmax\"]>", frac)
        label <- caption %>% strsplit("~") %>% getElement(1) %>% paste(collapse = " ")
        
        PM <- new("PMobj")
        PM@Name <- label
        PM@Caption <- caption
        PM@Stat <- xout/ref
        PM@Ref <- frac
        PM@Prob <- apply(PM@Stat[, yind, drop = FALSE] > PM@Ref, 2, mean)
        PM@Mean <- mean(PM@Stat[, yind, drop = FALSE])
        PM@MPs <- "Hist"
        PM
      })
    } else {
      PM <- RPC::make_PMobj(x = MSE, type = "SSB50%Rmax", frac = frac, year_range = dots$year_range)
    }
    
  } else if(type == "90%RS") {
    
    if(Hist) {
      PM <- local({
        yr <- seq(MSE@OM$CurrentYr[1] - MSE@nyears + 1, MSE@OM$CurrentYr[1])
        yind <- match(dots$year_range[1]:dots$year_range[2], yr)
        xout <- MSE@SSB_hist
        ref <- local({
          out <- RPC:::stock_recruit_int(MSE@Hist)
          RpS_90 <- apply(out$R/out$SSB, 1, quantile, probs = 0.9)
          R_90 <- apply(out$R, 1, quantile, probs = 0.9)
          R_90/RpS_90
        })
        caption <- paste0("Probability~SSB/SSB[\"90%R/S\"]>", frac)
        label <- caption %>% strsplit("~") %>% getElement(1) %>% paste(collapse = " ")
        
        PM <- new("PMobj")
        PM@Name <- label
        PM@Caption <- caption
        PM@Stat <- xout/ref
        PM@Ref <- frac
        PM@Prob <- apply(PM@Stat[, yind, drop = FALSE] > PM@Ref, 2, mean)
        PM@Mean <- mean(PM@Stat[, yind, drop = FALSE])
        PM@MPs <- "Hist"
        PM
      })
    } else {
      PM <- RPC::make_PMobj(x = MSE, type = "SSB90%R/S", frac = frac, year_range = dots$year_range)
    }
  } else if(type == "iSCAM_SSB0") {
    PM <- iSCAM_SSB0(MSE, frac = frac, year_range = dots$year_range, Hist = Hist)
  }
  return(PM)
}



#OM <- testOM
#OM@SRrel <- 2
#res <- runMSE(OM, MPs = c("FMSYref", "NFref"), extended = TRUE)
#SHPM(res, "HistSSB")
#SHPM(res, "50%Rmax")
#SHPM(res, "90%RS")
#SHPM(res, "SSBMSY")
#SHPM(res, "depletion", frac = 0.4)

FMSY_MP <- function(x, Data, reps = 1, frac = 1) {
  y <- max(Data@Year) - Data@LHYear+1
  nyears <- length(Data@Misc$FleetPars$Find[x,])
  FMSY <- Data@Misc$ReferencePoints$ByYear$FMSY[x,nyears+y]
  q <- Data@Misc$FleetPars$qs[x]
  qvar <- Data@Misc$FleetPars$qvar[x,y] # future only
  if (length(qvar)<1) qvar <- 1
  qinc <- Data@Misc$FleetPars$qinc[x] # future only
  qcur <- qvar * q*(1+qinc/100)^y # catchability this year
  
  HistE <- Data@OM$FinF[x] # Last historical fishing effort
  MSYE <- FMSY/qcur # effort for this year's FMSY
  
  Rec <- new('Rec')
  Rec@Effort <- MSYE/HistE * frac
  Rec
}

relF <- seq(0, 1.3, 0.1)

lapply(relF, function(x) {
  MP <- FMSY_MP
  formals(MP)$frac <- x
  MP_name <- paste0(x, "_FMSY")
  
  assign(MP_name, structure(MP, class = "MP"), envir = .GlobalEnv)
})

iSCAM_SSB0 <- function(MSE, frac, year_range, Hist = FALSE) {
  
  SSB0 <- local({
    nyears <- MSE@nyears
    maxage <- MSE@Hist@OM@maxage
    
    M <- apply(MSE@Hist@SampPars$Stock$M_ageArray[, , 1:nyears], 1:2, mean)
    Wt <- apply(MSE@Hist@SampPars$Stock$Wt_age[, , 1:nyears], 1:2, mean)
    Mat <- apply(MSE@Hist@SampPars$Stock$Mat_age[, , 1:nyears], 1:2, mean)
    Fec <- apply(MSE@Hist@SampPars$Stock$Fec_Age[, , 1:nyears], 1:2, mean)
    V <- rep(1, maxage + 1)
    R0 <- MSE@Hist@SampPars$Stock$R0
    h <- MSE@Hist@SampPars$Stock$hs
    SRrel <- MSE@Hist@OM@SRrel
    phi <- MSE@Hist@SampPars$Stock$SSBpR[, 1]
    
    sapply(1:MSE@nsim, function(x) {
      MSEtool:::MSYCalcs(logF = 0, M_at_Age = M[x, ], Wt_at_Age = Wt[x, ], Mat_at_Age = Mat[x, ], Fec_at_Age = Fec[x, ], V_at_Age = V, 
                         maxage = maxage, R0x = R0[x], SRrelx = SRrel, hx = h[x], SSBpR = phi[x], opt = 0L)["SB0"]
    })
  })
  caption <- paste0("Probability~SSB/~\"iSCAM\"~SSB[0]>", frac)
  label <- caption %>% strsplit("~") %>% getElement(1) %>% paste(collapse = " ")
  ref <- SSB0
  
  if(Hist) {
    yr <- seq(MSE@OM$CurrentYr[1] - MSE@nyears + 1, MSE@OM$CurrentYr[1])
    yind <- match(year_range[1]:year_range[2], yr)
    xout <- MSE@SSB_hist
    
    PM <- new("PMobj")
    PM@Name <- label
    PM@Caption <- caption
    PM@Stat <- xout/ref
    PM@Ref <- frac
    PM@Prob <- apply(PM@Stat[, yind, drop = FALSE] > PM@Ref, 2, mean)
    PM@Mean <- mean(PM@Stat[, yind, drop = FALSE])
    PM@MPs <- "Hist"
    
  } else {
    xout <- MSE@SSB
    
    CurrentYr <- MSE@OM$CurrentYr[1]
    
    nyp<-MSE@proyears
    nyh<-MSE@nyears
    pyind<-nyh+(1:nyp)
    py<-CurrentYr + 1:nyp
    
    yind <- match(year_range[1]:year_range[2], py)
    
    PM <- new("PMobj")
    PM@Name <- label
    PM@Caption <- caption
    PM@Stat <- xout/ref
    PM@Ref <- frac
    PM@Prob <- apply(PM@Stat[, , yind, drop = FALSE] > PM@Ref, c(1, 2), mean)
    PM@Mean <- apply(PM@Stat[, , yind, drop = FALSE] > PM@Ref, 2, mean)
    PM@MPs <- MSE@MPs
  }
  
  return(PM)
}

#B0 <- iSCAM_SSB0(MSE)
#t(SSB/B0) %>% matplot(type = 'l')
#abline(h = 0.3, lty = 3, col = 'red', lwd = 3)
