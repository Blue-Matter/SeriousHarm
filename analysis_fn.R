
get_ref_pt <- function(Hist, OM, Year_assess, 
                       Name = "",
                       Yr_SRP = Year[1],  # Year corresponding to ICES Blim
                       Yr_Fmed = Year[1], # First year for calculating median R/S
                       Yr_SP = NA,
                       Yr_Brecover,
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
  
  do_SP <- RPC::LRP_SP(Hist, "none")
  
  SP <- do_SP$Quantile[, "Median"]
  Removals <- do_SP$Removals[1, ]
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
                                relRfun = function(...) NULL,
                                SRRpars = list(),
                                maxage = OM@maxage,
                                plusgroup = 1)
      
      RPS <- Ref_search[3,]
      MSEtool:::LinInterp_cpp(RPS, F_search, xlev = median(RpS[Year >= Yr_Fmed]))
    })
  }
  
  
  if(missing(Yr_Brecover)) Yr_Brecover <- Year[which.min(SB)]
  if(!is.na(Yr_Brecover)) {
    Brecover <- SB[Year == Yr_Brecover]
  } else {
    Brecover <- NA
  }
  
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
                     SB_90 = SB90, SB_SP = B_SP, phi0 = phi0,
                     M = Hist@SampPars$Stock$Marray[1, 1:length(Year)],
                     Removals = Removals) %>% 
    mutate(Stock = Name)
  
  RPvars <- data.frame(RpS_med = median(RpS[Year >= Yr_Fmed]), RpS_90 = RpS_90, R_90 = R_90, 
                       Yr_Brecover = Yr_Brecover, Rmax50 = Rmax$Rmax50[1], 
                       SBRmax50 = SB_Rmax, Stock = Name, Yr_Fmed = Yr_Fmed, Yr_SP = Yr_SP)
  
  Assess <- local({
    a1 <- RPts[Year == Year_assess, ]
    cbind(data.frame(Year = Year_assess), a1$FM/a1[, c("FMSY", "Fmed")], a1["SPR"], a1["SPRdyn"], 
          a1$SB/a1[, c("SBMSY", "SB0", "SB0dyn", "SB_recover", "SB_Rmax", "SB_90", "SB_SP")])
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
    
    g
  })
  
  ggpubr::ggarrange(plotlist = out, legend = "bottom", common.legend = TRUE)
}



gghist <- function(x, title = TRUE, letter.legend, retro_only = TRUE, SBref = TRUE) {
  v <- x$RPts %>% 
    select(Year, SB, R, FM, M, Stock, FMSY) %>% 
    rename(Recruitment = R, F = FM) %>%
    reshape2::melt(id.vars = c("Year", "Stock")) %>%
    mutate(Mort = ifelse(variable == "M", "M", ifelse(variable == "FMSY", "FMSY", "F")), 
           variable = as.character(variable),
           linewidth = ifelse(variable == "SB", 1.5, 0.5))
  v$variable[v$variable %in% c("M", "F", "FMSY")] <- "Mortality~rate"
  v$variable <- factor(v$variable, levels = c("SB", "Mortality~rate", "Recruitment"))
  if(retro_only) v <- filter(v, Year <= x$Assess$Year)
  if(!title) v$Stock <- ""
  
  g <- ggplot(v, aes(Year, value))
  
  if (SBref) {
    SB_ref <- x$RPts %>%
      select(Year, SBMSY, SB0, SB0dyn, Stock) %>%
      reshape2::melt(id.var = c("Year", "Stock"), variable.name = "SBtype") %>%
      mutate(variable = factor("SB"))
    
    if(retro_only) {
      SB_ref <- filter(SB_ref, Year <= x$Assess$Year)
    }
    
    g <- g + 
      geom_line(data = SB_ref, aes(Year, value, group = SBtype), inherit.aes = FALSE) +
      geom_point(data = SB_ref %>% filter(Year %in% floor(seq(min(Year), max(Year), length.out = 10))), 
                 aes(Year, value, shape = SBtype), inherit.aes = FALSE) +
      scale_shape_manual("SB reference\npoint",
                         values = c("SB0" = 16, "SB0dyn" = 1, "SBMSY" = 8),
                         labels = c("SB0" = expression(SB[0~eq]), "SB0dyn" = expression(SB[0~dyn]), "SBMSY" = expression(SB[MSY])))
  }
  
  g <- g +
    geom_line(aes(linetype = Mort, linewidth = linewidth)) + 
    facet_grid(rows = vars(variable), 
               cols = vars(Stock), 
               labeller = label_parsed, 
               scales = "free", 
               switch = "y") +
    theme_bw() +
    expand_limits(y = 0) +
    theme(strip.background = element_rect(fill = NA, colour = NA), 
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "in"),
          strip.placement = "outside") +
    scale_linetype_manual("Mortality\nrate",
                          labels = c("F" = "F", "M" = "M", "FMSY" = expression(F[MSY])),
                          values = c("F" = 1, "M" = 2, "FMSY" = 3)) +
    scale_linewidth_identity()
  
  
  if(!retro_only) {
    g <- g + geom_vline(xintercept = x$Assess$Year, linetype = 2) 
  }
  
  if(!missing(letter.legend)) {
    lett <- data.frame(variable = factor(levels(v$variable), levels = levels(v$variable)),
                       Stock = v$Stock %>% unique(),
                       txt = letter.legend)
    
    g <- g + geom_text(data = lett, aes(label = txt), x = Inf, y = Inf, hjust = "inward", vjust = "inward")
  }
  
  g + guides(shape = guide_legend(order = 1), 
             linetype = guide_legend(order = 2))
  
}

plot_Hist <- function(...) {
  dots <- list(...)
  out <- lapply(dots, gghist)
  cowplot::plot_grid(plotlist = out)
}

ggSR <- function(x, year = TRUE, title = TRUE) {
  v <- x$RPts %>% select(Year, SB, R, Stock) 
  if(!title) v$Stock <- ""
  
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
}

plot_SR <- function(..., year = TRUE) {
  dots <- list(...)
  out <- lapply(dots, ggSR, year = year)
  cowplot::plot_grid(plotlist = out)
}

ggSR_LRP <- function(x, year = FALSE, 
                     title = TRUE, 
                     panel_labs = c("SRR", "SRP"),
                     size = 2, box.padding = 0.1, min.segment.length = 0.3,
                     letter.legend) {
  v <- x$RPts %>% select(Year, SB, R, Stock) 
  
  vout <- rbind(v, v) %>% mutate(Panel = panel_labs %>% rep(each = nrow(v))) %>%
    filter(Year >= as.numeric(x$RPvars["Yr_Fmed"])) %>%
    #filter(Panel == panel_labs[1] | Year >= x$RPvars["Yr_Fmed"] %>% as.numeric()) %>%
    mutate(fill = ifelse(Year < x$Assess$Year, "grey", 
                         ifelse(Year == x$Assess$Year, "black", "white")),
           Panel = factor(Panel, levels = panel_labs))
  
  SRpred <- x$SRpred %>% 
    mutate(Stock = x$RPvars$Stock, 
           Panel = factor(panel_labs[1], levels = panel_labs)) %>%
    filter(SB <= max(vout$SB))
  
  SRPvars <- as.data.frame(x$RPvars) %>% mutate(Panel = factor(panel_labs[2], levels = panel_labs))
  SRPvars$SB_90 <- unique(x$RPts$SB_90)
  
  SRRvars <- as.data.frame(x$RPvars) %>% mutate(Panel = factor(panel_labs[1], levels = panel_labs))
  
  #if(!title) vout$Stock <- SRpred$Stock <- SRPvars$Stock <- SRRvars$Stock <- ""
  
  g <- ggplot(vout, aes(SB, R)) + 
    geom_point(aes(fill = fill), shape = 21) + 
    geom_segment(data = SRPvars,                                        # 90 percentile R/S, and S
                 inherit.aes = FALSE, 
                 x = 0, y = 0, 
                 linetype = 3,
                 aes(xend = SB_90, yend = R_90)) +                      # 90 percentile R/S, and S
    geom_segment(data = SRPvars, 
                 inherit.aes = FALSE, 
                 x = 0, linetype = 3,
                 aes(xend = SB_90, y = R_90, yend = R_90)) +            # 90 percentile R/S, and S
    geom_vline(data = SRPvars,
               linetype = 2, 
               aes(xintercept = SB_90)) + 
    geom_point(data = SRPvars,                                          # 90 percentile R/S, and S
               inherit.aes = FALSE,    
               shape = 15, size = 3,
               aes(x = SB_90, y = R_90)) + 
    geom_segment(data = SRRvars,                                        # SSB at 50% Rmax
                 inherit.aes = FALSE, 
                 x = 0, linetype = 3,
                 aes(xend = SBRmax50, y = Rmax50, yend = Rmax50)) + 
    geom_vline(data = SRRvars,                                          # SSB at 50% Rmax
               linetype = 2,
               aes(xintercept = SBRmax50)) + 
    geom_line(data = SRpred) +                                          # Stock-recruit relationship
    geom_point(data = SRRvars, 
               inherit.aes = FALSE, 
               shape = 15, size = 3,
               aes(x = SBRmax50, y = Rmax50)) + 
    expand_limits(x = 0) + 
    #coord_cartesian(ylim = c(0, 1.1 * max(vout$R))) +
    theme_bw() +
    labs(y = "Recruitment") + 
    theme(panel.spacing = unit(0, "in"),
          strip.background = element_rect(fill = NA, colour = NA), 
          strip.placement = "outside") +
    scale_fill_identity()
  
  if(!missing(letter.legend)) {
    lett <- data.frame(Panel = factor(panel_labs, levels = panel_labs),
                       Stock = vout$Stock %>% unique(),
                       txt = letter.legend)
    
    g <- g + geom_text(data = lett, aes(label = txt), x = Inf, y = Inf, hjust = "inward", vjust = "inward")
  }
  
  
  if(year) {
    g <- g + 
      ggrepel::geom_text_repel(aes(label = Year), 
                               size = size,
                               box.padding = box.padding, 
                               min.segment.length = min.segment.length)
  }
  
  if(title) {
    g <- g +
      facet_grid(cols = vars(Stock), 
                 rows = vars(Panel), 
                 labeller = label_parsed)
  } else {
    g <- g + 
      facet_wrap(vars(Panel), labeller = label_parsed, ncol = 1, strip.position = "right")
  }
  
  g
}


plot_SR_LRP <- function(..., year = FALSE, size = 2, box.padding = 0.1, min.segment.length = 0.3) {
  dots <- list(...)
  
  out <- lapply(dots, ggSR_LRP, 
                year = year, size = size, 
                box.padding = box.padding, 
                min.segment.length = min.segment.length)
  
  cowplot::plot_grid(plotlist = out)
}


ggSP <- function(x, year = TRUE, SPB = TRUE, title = TRUE, letter.legend) {
  v <- x$SP %>% 
    filter(!is.na(SP)) %>% 
    mutate(SP_B = SP/B) %>% 
    rename(`Surplus~production~(SP)` = SP, `SP/B` = SP_B) %>% 
    reshape2::melt(id.vars = c("Year", "B")) %>%
    mutate(Stock = x$RPvars$Stock, 
           fill = ifelse(Year < x$Assess$Year, "grey", 
                         ifelse(Year == x$Assess$Year, "black", "white")))
  
  if(SPB) {
    v <- filter(v, variable == "SP/B")
  } else {
    v <- filter(v, variable == "Surplus~production~(SP)")
  }
  #if(!title) v$Stock <- ""
  
  g <- ggplot(v, aes(B, value)) + 
    geom_path() + 
    geom_point(aes(fill = fill), shape = 21) +
    geom_hline(yintercept = 0, linetype = 3) + 
    #facet_grid(vars(variable), vars(Stock), labeller = label_parsed, scales = "free_y", switch = "y") + 
    expand_limits(x = 0, y = 0) + 
    theme_bw() +
    theme(strip.background = element_rect(fill = NA, colour = NA), 
          strip.placement = "outside",
          axis.title.y = element_blank(),
          panel.spacing.y = unit(0, "in")) + 
    labs(x = "Total Biomass (B)") +
    scale_fill_identity()
  
  if(!is.na(x$RPvars$Yr_SP)) {
    B_SP <- x$SP$B[x$SP$Year == x$RPvars$Yr_SP]
    g <- g + geom_vline(xintercept = B_SP, linetype = 2)
  }
  if(!missing(letter.legend)) {
    lett <- data.frame(variable = factor(levels(v$variable), levels = levels(v$variable)),
                       Stock = v$Stock %>% unique(),
                       txt = letter.legend)
    if(SPB) { # May not need
      lett <- filter(lett, variable == "SP/B")
    } else {
      lett <- filter(lett, variable == "Surplus~production~(SP)")
    }
    g <- g + geom_text(data = lett, aes(label = txt), x = Inf, y = Inf, hjust = "inward", vjust = "inward")
    
  }
  
  if(year) g <- g + ggrepel::geom_text_repel(aes(label = Year), max.overlaps = 5)
  
  if(title) {
    g <- g + 
      facet_grid(vars(variable), vars(Stock), labeller = label_parsed, scales = "free_y", switch = "y")
  } else {
    g <- g +
      facet_wrap(vars(variable), labeller = label_parsed, scales = "free_y", strip.position = "left")
  }
  
  g
}

plot_SP <- function(..., year = TRUE, SPB = TRUE) {
  dots <- list(...)
  
  out <- lapply(dots, ggSP, year = year, SPB = SPB)
  cowplot::plot_grid(plotlist = out)
}


SRR <- function(x, R, S, rel = c("hockey_stick", "Ricker", "BH", "sBH"), 
                opt = TRUE, figure = TRUE, delta) {
  
  rel <- match.arg(rel)
  
  if(rel == "hockey_stick") {
    Scp <- x[1] * max(S)
  } else {
    beta <- exp(x[1])/mean(S)
  }
  
  if(rel == "sBH" && missing(delta)) {
    delta <- x[2]
  }
  
  f_S <- switch(rel,
                "hockey_stick" = ifelse(S <= Scp, S, Scp),
                "Ricker" = S * exp(-beta * S),
                "BH" = S / (1 + beta * S),
                "sBH" = S ^ delta / (1 + (beta * S)^delta) # K = 1/beta
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
    
    if(rel == "sBH") {
      out$delta <- delta
    }
    
    if(figure) {
      plot(S, R, xlab = "Spawners", ylab = "Recruitment", xlim = c(0, 1.1 * max(S)), ylim = c(0, 1.1 * max(R)))
      
      Splot <- seq(0, 1.1 * max(S), length.out = 100)
      Rplot <- alpha * switch(rel,
                              "hockey_stick" = ifelse(Splot <= Scp, Splot, Scp),
                              "Ricker" = Splot * exp(-beta * Splot),
                              "BH" = Splot / (1 + beta * Splot),
                              "sBH" = Splot^delta / (1 + (beta * Splot)^delta)
      )
      lines(Splot, Rplot, col = 'red')
    }
    
    return(out)
  }
}

calc_phi0 <- function(M, Wt, Mat, Fec) {
  1/MSEtool:::Ref_int_cpp(1e-8, M, Wt, Mat, Fec, rep(0, length(M)),
                          relRfun = function(...) NULL, 
                          SRRpars = list(),
                          length(M) - 1)[3, 1]
}

plot_ts <- function(LRP) {
  
  var_order <- c("F/F[MSY]", "SPR[eq]", "SPR[dyn]", "SB/SB[MSY]", "SB/SB[0~eq]", "SB/SB[0~dyn]",
                 "SB/SB[0.5~Rmax]", "SB/SB[low~SP]", "SB/SB[recover]", "SB/SB[0.9~R/S]")
  v <- LRP$RPts %>% group_by(Year) %>%
    mutate(`F/F[MSY]` = FM/FMSY, 
           `SB/SB[MSY]` = SB/SBMSY, `SB/SB[0~eq]` = SB/SB0, `SB/SB[0~dyn]` = SB/SB0dyn,
           `SPR[eq]` = SPR, `SPR[dyn]` = SPRdyn, `SB/SB[recover]` = SB/SB_recover,
           `SB/SB[0.5~Rmax]` =  SB/SB_Rmax,
           `SB/SB[0.9~R/S]` = SB/SB_90,
           `SB/SB[low~SP]` = SB/SB_SP,
           .keep = "none"
    ) %>%
    #left_join(LRP$SP, by = "Year") %>% 
    #mutate(`SP/B` = SP/B) %>%
    #select(!B) %>%
    reshape2::melt(id.var = "Year") %>%
    mutate(variable = factor(variable, levels = var_order)) %>%
    filter(!is.na(value))
  
  ref_vars <- var_order[var_order %in% unique(v$variable)]
  ref_vars <- factor(ref_vars, levels = ref_vars)
  
  ref <- data.frame(variable = ref_vars) %>% 
    mutate(name = paste0("(", letters[1:nrow(.)], ")"))
  
  g <- ggplot(v, aes(Year, value)) + 
    geom_path() + 
    geom_text(data = ref, aes(label = name), x = -Inf, y = Inf, hjust = "inward", vjust = "inward") + 
    geom_vline(xintercept = LRP$Assess$Year, linetype = 2) + 
    facet_wrap(vars(variable), 
               labeller = label_parsed, 
               ncol = 3,
               nrow = 4,
               scales = "free_y", 
               strip.position = "left") + 
    expand_limits(y = 0) + 
    theme_bw() +
    theme(strip.background = element_rect(fill = NA, colour = NA), 
          strip.placement = "outside",
          axis.title.y = element_blank(),
          panel.spacing.y = unit(0, "in")) +
    ggtitle(parse(text = LRP$RPts$Stock %>% unique()))
  
  g
}

#plot_summary_ts <- function(LRP, vars = "FMSY") {
#  vars_all <- data.frame(vars = c("FMSY", "SPR", "SPRdyn", "SBMSY", "SB0", "SB0dyn", "SB_recover", "SB_Rmax", "SB_90", "SB_SP"),
#                         vtitle = c("F/F[MSY]", "SPR[eq]", "SPR[dyn]", "SB/SB[MSY]", "SB/SB[0~eq]", "SB/SB[0~dyn]", 
#                                    "SB/SB[recover]", "SB/SB[0.5~Rmax]", "SB/SB[0.9~R/S]", "SB/SB[low~SP]"),
#                         type = c("FM", "SPR", "SPR", "SB", "SB", "SB", "SB", "SB", "SB", "SB"))
#  
#  df <- lapply(vars, function(x) {
#    ind <- match(x, vars_all$vars)
#    
#    if(vars_all$type[ind] == "SPR") {
#      data.frame(Year = LRP$RPts$Year, 
#                 value = LRP$RPts[, vars_all$vars[ind]],
#                 vars = x, 
#                 type = vars_all$type[ind], 
#                 name = vars_all$vtitle[ind])
#    } else {
#      data.frame(Year = LRP$RPts$Year, 
#                 value = LRP$RPts[, vars_all$type[ind]] / LRP$RPts[, x],
#                 vars = x, 
#                 type = vars_all$type[ind], 
#                 name = vars_all$vtitle[ind])
#    }
#  }) %>% bind_rows()
#  
#  g <- ggplot(df, aes(Year, value)) + 
#    geom_line(aes(linetype = name)) + 
#    theme_bw() + 
#    expand_limits(y = 0)
#  
#  browser()
#  
#  if("SPR" %in% vars) {
#    g <- g + 
#      labs(y = "SPR") +
#      scale_linetype_manual(name = "SPR", values = c("SPR[eq]" = 1, "SPR[dyn]" = 2),
#                            labels = c("SPR[eq]" = "Equilibrium", "SPR[dyn]" = "Dynamic"))
#  } else { # SB
#    g <- g + 
#      labs(y = "Biomass LRP ratio") +
#      scale_linetype_manual(name = "LRP", values = label_value, 
#                            labels = label_parsed)
#  }
#  g
#}
#
#plot_summary_ts(LRP, c("SPR", "SPRdyn"))
#
#
#plot_summary_ts(LRP, c("SBMSY", "SB0"))
