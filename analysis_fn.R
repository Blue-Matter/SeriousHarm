
get_ref_pt <- function(Hist, OM, Year_assess, 
                       Name = "",
                       Yr_SRP = Year[1],  # Year corresponding to ICES Blim
                       Yr_Fmed = Year[1], # First year for calculating median R/S
                       Yr_SP = NA,
                       figure = FALSE, 
                       BSRP) {
  Year <- OM@CurrentYr - OM@nyears:1 + 1
  
  SB <- Hist@TSdata$SBiomass[1, , ] %>% rowSums()
  R <- Hist@AtAge$Number[1, 1, , ] %>% rowSums()
  RpS <- R/SB
  
  phi0 <- sapply(1:OM@nyears, function(y) {
    calc_phi0(M = Hist@SampPars$Stock$M_ageArray[1, , y], 
              Wt = Hist@SampPars$Stock$Wt_age[1, , y], 
              Mat = Hist@SampPars$Stock$Mat_age[1, , y], 
              Fec = Hist@SampPars$Stock$Fec_Age[1, , y])
  })
  
  
  SP <- RPC:::hist_BvsSP(Hist, FALSE)[, "Median"]
  B <- Hist@TSdata$Biomass[1, , ] %>% rowSums()
  
  SBMSY <- Hist@Ref$ByYear$SSBMSY[1, 1:OM@nyears]
  SB0 <- Hist@Ref$ByYear$SSB0[1, 1:OM@nyears]
  SB0dyn <- Hist@Ref$Dynamic_Unfished$SSB0[1, 1:OM@nyears]
  
  FM <- Hist@AtAge$F.Mortality[1, , , ] %>% apply(2, max)
  FMSY <- Hist@Ref$ByYear$FMSY[1, 1:OM@nyears]
  SPR <- Hist@TSdata$SPR$Equilibrium[1, ]
  SPRdyn <- Hist@TSdata$SPR$Dynamic[1, ]
  
  Fmed <- Hist@Ref$ByYear$Fmed[1, 1:OM@nyears]
  
  if(Yr_Fmed == Year[1]) {
    Fmed <- Hist@Ref$ByYear$Fmed[1, 1:OM@nyears]
  } else {
    Fmed <- sapply(1:OM@nyears, function(y) {
      boundsF <- c(1E-3, 3)
      F_search <- exp(seq(log(min(boundsF)), log(max(boundsF)), length.out = 50))
      
      Ref_search <- MSEtool:::Ref_int_cpp(F_search, M_at_Age = Hist@SampPars$Stock$M_ageArray[1, , y],
                                Wt_at_Age = Hist@SampPars$Stock$Wt_age[1, , y], 
                                Mat_at_Age = Hist@SampPars$Stock$Mat_age[1, , y],
                                Fec_at_Age = Hist@SampPars$Stock$Fec_Age[1, , y],
                                V_at_Age = Hist@SampPars$Fleet$V_real[1, , y],
                                maxage = OM@maxage,
                                plusgroup = 1)
      
      RPS <- Ref_search[3,]
      MSEtool:::LinInterp_cpp(RPS, F_search, xlev = median(RpS[Year >= Yr_Fmed]))
    })
  }
  
  
  Yr_Brecover <- Year[which.min(SB)]
  Brecover <- min(SB)
  
  if(missing(BSRP)) {
    if(is.character(Yr_SRP)) {
      if(Yr_SRP == "hockey_stick") {
        #browser()
        #x <- seq(0, 1, 0.01)
        #prof <- sapply(x, SRR, R = R[Year >= Yr_Fmed], S = SB[Year >= Yr_Fmed])
        #plot(x, prof)
        #vars <- SRR(x[which.min(prof)], R = R[Year >= Yr_Fmed], S = SB[Year >= Yr_Fmed], opt = FALSE)
        #vars$Scp
        
        opt <- optimize(SRR, interval = c(0, 1), R = R[Year >= Yr_Fmed], S = SB[Year >= Yr_Fmed])
        vars <- SRR(opt$minimum, R = R[Year >= Yr_Fmed], S = SB[Year >= Yr_Fmed], opt = FALSE)
        BSRP <- vars$Scp
      } else {
        stop("Invalid argument for Yr_SRP")
      }
    } else if(is.numeric(Yr_SRP) && !is.na(Yr_SRP)) {
      BSRP <- SB[which(Year == Yr_SRP)]
    } else {
      BSRP <- NA_real_
    }
  }
  
  
  if(!is.na(Yr_SP)) {
    B_SP <- SB[Year == Yr_SP]
  } else {
    B_SP <- NA
  }
  
  Rmax <- RPC:::calculate_SSB50(Hist)
  SB_Rmax <- Rmax$SSB50[1]
  
  SRpred <- RPC:::stock_recruit_int(Hist)
  
  RpS_90 <- quantile(RpS[Year >= Yr_Fmed], probs = 0.9)
  R_90 <- quantile(R[Year >= Yr_Fmed], probs = 0.9)
  SB90 <- R_90/RpS_90
  
  RPts <- data.frame(Year = Year, FM = FM, FMSY = FMSY, Fmed = Fmed, SPR = SPR, SPRdyn = SPRdyn, R = R, 
                     SB = SB, SBMSY = SBMSY, SB0 = SB0, SB0dyn = SB0dyn, SB_SRP = BSRP, SB_recover = Brecover, SB_Rmax = SB_Rmax, 
                     SB_90 = SB90, SB_SP = B_SP, phi0 = phi0) %>% mutate(Stock = Name)
  
  RPvars <- data.frame(RpS_med = median(RpS[Year >= Yr_Fmed]), RpS_90 = RpS_90, R_90 = R_90, 
                       Yr_Brecover = Yr_Brecover, Rmax50 = Rmax$Rmax50[1], 
                       SBRmax50 = SB_Rmax, Stock = Name, Yr_Fmed = Yr_Fmed, Yr_SP = Yr_SP)
  
  Assess <- local({
    a1 <- RPts[Year == Year_assess, ]
    cbind(data.frame(Year = Year_assess), a1$FM/a1[, 3:4], a1["SPR"], a1["SPRdyn"], a1$SB/a1[, c(9:16)])
  })
  print(Assess)
  print(RPvars)
  invisible(list(Assess = Assess, RPts = RPts, RPvars = RPvars, 
                 SP = data.frame(Year = Year, B = B, SP = c(SP, NA)),
                 SRpred = data.frame(SB = SRpred$predSSB, R = SRpred$predR[1, ])))
}


plot_F <- function(...) {
  dots <- list(...)
  
  out <- lapply(1:length(dots), function(i) {
    v <- dots[[i]]$RPts %>% 
      mutate(F_FMSY = FM/FMSY, F_Fmed = FM/Fmed) %>% 
      select(Year, SPR, FM, F_FMSY, F_Fmed, Stock) %>% 
      rename(F = FM, `F/F[MSY]` = F_FMSY, `F/F[med]` = F_Fmed) %>%
      reshape2::melt(id.vars = c("Year", "Stock")) %>% mutate(Type = "Equilibrium")
    
    v <- rbind(v,
               data.frame(Year = unique(v$Year), Stock = unique(v$Stock),
                          variable = "SPR", value = dots[[i]]$RPts$SPRdyn, Type = "Dynamic"))
    
    yint <- data.frame(yy = c(0, 1, 0, 0), variable = c("F", "SPR", "F/F[MSY]", "F/F[med]")) %>%
      mutate(variable = factor(variable, levels = c("F", "SPR", "F/F[MSY]", "F/F[med]")))
    
    g <- ggplot(v, aes(Year, value)) + 
      geom_line(aes(linetype = Type)) + 
      facet_grid(rows = vars(variable), 
                 labeller = label_parsed,
                 cols = vars(Stock), 
                 scales = "free", switch = "y") + 
      geom_blank(data = yint, inherit.aes = FALSE, aes(y = yy)) + 
      geom_vline(xintercept = dots[[i]]$Assess$Year, linetype = 3) +
      expand_limits(y = 0) + 
      theme_bw() +
      scale_linetype_manual(name = "SPR Type", values = c("Equilibrium" = 1, "Dynamic" = 2)) +
      theme(strip.background = element_rect(fill = NA, colour = NA), 
            axis.title.y = element_blank(),
            panel.spacing = unit(0, "in"),
            strip.placement = "outside")
    #if(i > 1) g <- g + theme(strip.text.y = element_blank())
    g
  })
  
  ggpubr::ggarrange(plotlist = out, legend = "bottom", common.legend = TRUE)
}


plot_B <- function(...) {
  dots <- list(...)
  
  out <- lapply(1:length(dots), function(i) {
    v <- dots[[i]]$RPts %>% 
      mutate(SB_SBMSY = SB/SBMSY, SB_SB0 = SB/SB0) %>% 
      select(Year, SB, SB_SBMSY, SB_SB0, Stock) %>% 
      rename(`SB/SB[MSY]` = SB_SBMSY, `SB/SB[0]` = SB_SB0) %>%
      reshape2::melt(id.vars = c("Year", "Stock")) %>%
      mutate(Type = "Equilibrium")
    
    v <- rbind(v,
               data.frame(Year = unique(v$Year), Stock = unique(v$Stock),
                          variable = "SB/SB[0]", value = dots[[i]]$RPts$SB/dots[[i]]$RPts$SB0dyn, Type = "Dynamic"))
    
    g <- ggplot(v, aes(Year, value)) + 
      geom_line(aes(linetype = Type)) + 
      facet_grid(rows = vars(variable), 
                 labeller = label_parsed,
                 cols = vars(Stock), 
                 scales = "free", switch = "y") + 
      expand_limits(y = 0) + 
      geom_vline(xintercept = dots[[i]]$Assess$Year, linetype = 3) +
      theme_bw() +
      scale_linetype_manual(name = parse(text = "SB[0]~Type"), values = c("Equilibrium" = 1, "Dynamic" = 2)) +
      theme(strip.background = element_rect(fill = NA, colour = NA), 
            axis.title.y = element_blank(),
            panel.spacing = unit(0, "in"),
            strip.placement = "outside")
    #if(i > 1) g <- g + theme(strip.text.y = element_blank())
    g
  })
  
  ggpubr::ggarrange(plotlist = out, legend = "bottom", common.legend = TRUE)
}





plot_Hist <- function(...) {
  dots <- list(...)
  
  out <- lapply(1:length(dots), function(i) {
    v <- dots[[i]]$RPts %>% 
      select(Year, SB, R, FM, Stock) %>% 
      rename(Recruitment = R, F = FM) %>%
      reshape2::melt(id.vars = c("Year", "Stock")) 
    
    ggplot(v, aes(Year, value)) + 
      geom_line() + 
      facet_grid(rows = vars(variable), cols = vars(Stock), 
                 labeller = label_parsed, scales = "free", switch = "y") + 
      expand_limits(y = 0) + 
      geom_vline(xintercept = dots[[i]]$Assess$Year, linetype = 2) + 
      theme_bw() +
      theme(strip.background = element_rect(fill = NA, colour = NA), 
            axis.title.y = element_blank(),
            panel.spacing = unit(0, "in"),
            strip.placement = "outside")
  })
  
  cowplot::plot_grid(plotlist = out)
}


plot_SR <- function(..., year = TRUE) {
  dots <- list(...)
  
  out <- lapply(1:length(dots), function(i) {
    v <- dots[[i]]$RPts %>% 
      select(Year, SB, R, Stock) 
    
    g <- ggplot(v, aes(SB, R)) + 
      geom_point(shape = 21, fill = "grey") + 
      facet_grid(cols = vars(Stock), labeller = label_parsed, scales = "free", switch = "y") + 
      expand_limits(x = 0, y = 0) + 
      theme_bw() +
      labs(y = "Recruitment") + 
      theme(strip.background = element_rect(fill = NA, colour = NA), 
            strip.placement = "outside",
            panel.spacing.y = unit(0, "in"))
    
    if(year) {
      g <- g + ggrepel::geom_text_repel(aes(label = Year))
    }
    g
  })
  
  cowplot::plot_grid(plotlist = out)
}



plot_SR_LRP <- function(..., year = FALSE) {
  dots <- list(...)
  
  out <- lapply(1:length(dots), function(i) {
    
    v <- dots[[i]]$RPts %>% 
      select(Year, SB, R, Stock) 
    
    vout <- rbind(v, v) %>% mutate(Panel = c("SRP", "SRR") %>% rep(each = nrow(v))) %>%
      filter(Panel == "SRR" | Year >= dots[[i]]$RPvars["Yr_Fmed"] %>% as.numeric())
    
    SRpred <- dots[[i]]$SRpred %>% mutate(Stock = dots[[i]]$RPvars$Stock, Panel = "SRR")
    
    SRPvars <- as.data.frame(dots[[i]]$RPvars) %>% mutate(Panel = "SRP")
    SRPvars$SB_90 <- unique(dots[[i]]$RPts$SB_90)
    
    SRRvars <- as.data.frame(dots[[i]]$RPvars) %>% mutate(Panel = "SRR")
    
    g <- ggplot(vout, aes(SB, R)) + 
      geom_point(shape = 21, fill = "grey") + 
      geom_abline(data = SRPvars, 
                  aes(intercept = 0, slope = RpS_med)) + 
      geom_segment(data = SRPvars, inherit.aes = FALSE, x = 0, y = 0, linetype = 2,
                   aes(xend = SB_90, yend = R_90)) + 
      geom_segment(data = SRPvars, inherit.aes = FALSE, x = 0, linetype = 2,
                   aes(xend = SB_90, y = R_90, yend = R_90)) + 
      geom_point(data = SRPvars, inherit.aes = FALSE, shape = 15, size = 3, #col = "red",
                 aes(x = SB_90, y = R_90)) + 
      geom_segment(data = SRRvars, inherit.aes = FALSE, x = 0, linetype = 3,
                   aes(xend = SBRmax50, y = Rmax50, yend = Rmax50)) + 
      geom_segment(data = SRRvars, inherit.aes = FALSE, y = 0, linetype = 3,
                   aes(x = SBRmax50, xend = SBRmax50, yend = Rmax50)) + 
      geom_line(data = SRpred) + 
      geom_point(data = SRRvars, inherit.aes = FALSE, shape = 15, size = 3, #col = "red",
                   aes(x = SBRmax50, y = Rmax50)) + 
      facet_grid(cols = vars(Stock), rows = vars(Panel), labeller = label_parsed, scales = "free_x") + 
      expand_limits(x = 0, y = 0) + 
      theme_bw() +
      labs(y = "Recruitment") + 
      theme(strip.background = element_rect(fill = NA, colour = NA), 
            strip.placement = "outside")
            #panel.spacing.y = unit(0, "in"))
    
    if(year) g <- g + ggrepel::geom_text_repel(aes(label = Year))
    g
  })
  
  cowplot::plot_grid(plotlist = out)
}


plot_SP <- function(..., year = TRUE) {
  dots <- list(...)
  
  out <- lapply(1:length(dots), function(i) {
    
    v <- dots[[i]]$SP %>% filter(!is.na(SP)) %>% mutate(SP_B = SP/B) %>% 
      rename(`Surplus~production~(SP)` = SP, `SP/B` = SP_B) %>% reshape2::melt(id.vars = c("Year", "B")) %>%
      mutate(Stock = dots[[i]]$RPvars$Stock)
    
    g <- ggplot(v, aes(B, value)) + 
      geom_path() + 
      geom_point(shape = 21, fill = "grey") +
      geom_hline(yintercept = 0, linetype = 2) + 
      facet_grid(vars(variable), vars(Stock), labeller = label_parsed, scales = "free_y", switch = "y") + 
      expand_limits(x = 0, y = 0) + 
      theme_bw() +
      theme(strip.background = element_rect(fill = NA, colour = NA), 
            strip.placement = "outside",
            axis.title.y = element_blank(),
            panel.spacing.y = unit(0, "in")) + 
      labs(x = "Biomass (B)")
    if(year) g <- g + ggrepel::geom_text_repel(aes(label = Year))
    g
  })
  
  cowplot::plot_grid(plotlist = out)
}


SRR <- function(x, R, S, rel = c("hockey_stick", "Ricker", "BH"), opt = TRUE, figure = TRUE) {
  
  rel <- match.arg(rel)
  
  if(rel == "hockey_stick") {
    Scp <- x[1] * max(S)
  } else {
    beta <- exp(x[1])/mean(S)
  }
  
  f_S <- switch(rel,
                "hockey_stick" = ifelse(S <= Scp, S, Scp),
                "Ricker" = S * exp(-beta * S),
                "BH" = S / (1 + beta * S)
  )
  
  alpha <- mean(log(R/f_S)) %>% exp()
  Rpred <- alpha * f_S
  
  sigmaR <- mean(log(R/Rpred) * log(R/Rpred)) %>% sqrt()
  
  if(opt) {
    
    nll <- dnorm(log(R), log(Rpred), sigmaR, log = TRUE) %>% sum()
    return(-1 * nll)
    
  } else {
    
    out <- list(alpha = alpha, Rpred = Rpred, sigmaR = sigmaR)
    if(rel == "hockey_stick") {
      out$Scp <- Scp
      out$Rmax <- alpha * Scp
    } else {
      out$beta <- beta
    }
    
    if(figure) {
      plot(S, R, xlab = "Spawners", ylab = "Recruitment", xlim = c(0, 1.1 * max(S)), ylim = c(0, 1.1 * max(R)))
      
      Splot <- seq(0, 1.1 * max(S), length.out = 100)
      Rplot <- alpha * switch(rel,
                      "hockey_stick" = ifelse(Splot <= Scp, Splot, Scp),
                      "Ricker" = Splot * exp(-beta * Splot),
                      "BH" = Splot / (1 + beta * Splot)
      )
      lines(Splot, Rplot, col = 'red')
    }
    
    return(out)
  }
}

calc_phi0 <- function(M, Wt, Mat, Fec) {
  1/MSEtool:::Ref_int_cpp(1e-8, M, Wt, Mat, Fec, rep(0, length(M)), length(M) - 1)[3, 1]
}

