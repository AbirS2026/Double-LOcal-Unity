# ============================================================
# Camera-ready code for Figure 3 in the main paper
#
# Figure 3: Estimated autoregressive coefficients and
#           local-deviation parameters across asset classes
#
# Panels:
#   (a) Commodity markets, 2022--2026
#       Output: figures/commodity_delta_gamma_two_regimes.pdf
#
#   (b) Global real-estate markets, 2020--2026
#       Output: figures/global_real_estate_delta_gamma_two_regimes.pdf
#
# The figure reports \hat{rho}_n and \hat{gamma}_n estimates with
# asymptotic 95% confidence intervals. The dashed horizontal line
# marks unity.
# ============================================================
# ============================================================


library(quantmod)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# ------------------------------------------------------------------------------
# Global options
# ------------------------------------------------------------------------------

SAVE_FIGURES <- TRUE
SAVE_TABLES  <- TRUE

OUTPUT_DIR <- "figures"
TABLE_DIR  <- "tables"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(TABLE_DIR)) dir.create(TABLE_DIR, recursive = TRUE)

# ------------------------------------------------------------------------------
# Asset lists
# ------------------------------------------------------------------------------

COMMODITY_ASSETS <- data.frame(
  Asset = c(
    "Lithium ETF",
    "Crude Oil",
    "Rare Earth ETF",
    "Copper",
    "Gold",
    "Uranium ETF",
    "Silver",
    "Natural Gas"
  ),
  Ticker = c(
    "LIT",
    "CL=F",
    "REMX",
    "HG=F",
    "GC=F",
    "URA",
    "SI=F",
    "NG=F"
  ),
  stringsAsFactors = FALSE
)

GLOBAL_RE_ASSETS <- data.frame(
  Market = c(
    "United States",
    "Canada",
    "Australia",
    "United Kingdom",
    "Japan",
    "India",
    "China"
  ),
  Ticker = c(
    "VNQ",
    "XRE.TO",
    "VAP.AX",
    "IUKP.L",
    "1343.T",
    "DLF.NS",
    "1109.HK"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Data download
# ------------------------------------------------------------------------------

get_yahoo <- function(ticker, label, group_name, start_date, end_date, min_obs = 100) {
  x <- tryCatch(
    quantmod::getSymbols(
      ticker,
      src = "yahoo",
      from = as.Date(start_date),
      to = as.Date(end_date),
      auto.assign = FALSE,
      warnings = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(x)) {
    warning(paste("Could not fetch", label, ticker))
    return(NULL)
  }
  
  px <- tryCatch(quantmod::Ad(x), error = function(e) quantmod::Cl(x))
  
  out <- data.frame(
    Date = as.Date(index(px)),
    Ticker = ticker,
    Value = as.numeric(px),
    stringsAsFactors = FALSE
  )
  
  out[[group_name]] <- label
  
  out <- out |>
    na.omit() |>
    select(Date, all_of(group_name), Ticker, Value)
  
  if (nrow(out) < min_obs) {
    warning(paste("Too few observations for", label, ticker))
    return(NULL)
  }
  
  out
}

make_asset_data <- function(asset_df, group_name, start_date, end_date) {
  out <- bind_rows(lapply(seq_len(nrow(asset_df)), function(i) {
    get_yahoo(
      ticker = asset_df$Ticker[i],
      label = asset_df[[group_name]][i],
      group_name = group_name,
      start_date = start_date,
      end_date = end_date
    )
  }))
  
  out |>
    arrange(.data[[group_name]], Date)
}

# ------------------------------------------------------------------------------
# Delta-gamma estimation
# ------------------------------------------------------------------------------

compute_delta_gamma_stats <- function(X, lambda = 0.05, normalize_start = TRUE) {
  X <- as.numeric(X)
  X <- X[is.finite(X)]
  
  if (length(X) < 3) {
    return(c(
      T = NA_real_, delta_hat = NA_real_, gamma_hat = NA_real_,
      delta_lower = NA_real_, delta_upper = NA_real_,
      gamma_lower = NA_real_, gamma_upper = NA_real_,
      regime_code = NA_real_
    ))
  }
  
  if (normalize_start) X <- X / X[1]
  
  y <- X[-1]
  xlag <- X[-length(X)]
  Tn <- length(y)
  
  delta_hat <- sum(y * xlag) / sum(xlag^2)
  diff_from_one <- abs(delta_hat - 1)
  
  gamma_hat <- if (diff_from_one <= .Machine$double.eps) {
    NA_real_
  } else {
    -log(diff_from_one) / log(Tn)
  }
  
  cv_normal <- qnorm(1 - lambda / 2)
  cv_cauchy <- qcauchy(1 - lambda / 2)
  
  if (is.na(gamma_hat)) {
    return(c(
      T = Tn,
      delta_hat = delta_hat,
      gamma_hat = NA_real_,
      delta_lower = NA_real_,
      delta_upper = NA_real_,
      gamma_lower = NA_real_,
      gamma_upper = NA_real_,
      regime_code = NA_real_
    ))
  }
  
  if (delta_hat < 1) {
    regime_code <- 0
    delta_hw <- cv_normal * sqrt(2) / Tn^((1 + gamma_hat) / 2)
    gamma_hw <- cv_normal * sqrt(2) / (Tn^((1 - gamma_hat) / 2) * log(Tn))
  } else {
    regime_code <- 1
    delta_hw <- cv_cauchy * 2 / ((Tn^gamma_hat) * (delta_hat^Tn))
    gamma_hw <- cv_cauchy * 2 / (((1 + 1 / Tn^gamma_hat)^Tn) * log(Tn))
  }
  
  c(
    T = Tn,
    delta_hat = delta_hat,
    gamma_hat = gamma_hat,
    delta_lower = delta_hat - delta_hw,
    delta_upper = delta_hat + delta_hw,
    gamma_lower = gamma_hat - gamma_hw,
    gamma_upper = gamma_hat + gamma_hw,
    regime_code = regime_code
  )
}

# ------------------------------------------------------------------------------
# Camera-ready plot function
# ------------------------------------------------------------------------------

make_delta_gamma_plot <- function(
    asset_df,
    group_name,
    start_date,
    end_date,
    lambda = 0.05,
    table_file,
    figure_file_base,
    save_figure = SAVE_FIGURES,
    save_tables = SAVE_TABLES
) {
  
  # -----------------------------
  # Download and estimate
  # -----------------------------
  
  price_df <- make_asset_data(
    asset_df = asset_df,
    group_name = group_name,
    start_date = start_date,
    end_date = end_date
  )
  
  est_df <- price_df |>
    group_by(.data[[group_name]], Ticker) |>
    summarise(
      stats = list(compute_delta_gamma_stats(Value, lambda = lambda)),
      .groups = "drop"
    ) |>
    mutate(
      T = sapply(stats, `[[`, "T"),
      Delta_Hat = sapply(stats, `[[`, "delta_hat"),
      Gamma_Hat = sapply(stats, `[[`, "gamma_hat"),
      Delta_CI_Lower = sapply(stats, `[[`, "delta_lower"),
      Delta_CI_Upper = sapply(stats, `[[`, "delta_upper"),
      Gamma_CI_Lower = sapply(stats, `[[`, "gamma_lower"),
      Gamma_CI_Upper = sapply(stats, `[[`, "gamma_upper"),
      Regime_Code = sapply(stats, `[[`, "regime_code"),
      Regime = ifelse(
        Regime_Code == 1,
        "Mildly explosive",
        "Mildly integrated"
      )
    ) |>
    select(-stats)
  
  est_df[[group_name]] <- factor(
    est_df[[group_name]],
    levels = asset_df[[group_name]]
  )
  
  # -----------------------------
  # Plot data
  # -----------------------------
  
  plot_df <- bind_rows(
    est_df |>
      transmute(
        Item = .data[[group_name]],
        Ticker,
        Regime,
        Metric = "Estimated Autoregressive Coefficient",
        Estimate = Delta_Hat,
        Lower = Delta_CI_Lower,
        Upper = Delta_CI_Upper,
        Reference = 1
      ),
    est_df |>
      transmute(
        Item = .data[[group_name]],
        Ticker,
        Regime,
        Metric = "Estimated Growth Rate",
        Estimate = Gamma_Hat,
        Lower = Gamma_CI_Lower,
        Upper = Gamma_CI_Upper,
        Reference = 1
      )
  )
  
  plot_df$Item <- factor(plot_df$Item, levels = asset_df[[group_name]])
  
  plot_df$Metric <- factor(
    plot_df$Metric,
    levels = c(
      "Estimated Autoregressive Coefficient",
      "Estimated Growth Rate"
    )
  )
  
  plot_df$Regime <- factor(
    plot_df$Regime,
    levels = c("Mildly explosive", "Mildly integrated")
  )
  
  # -----------------------------
  # Colors and shapes
  # -----------------------------
  
  regime_cols <- c(
    "Mildly explosive"  = "#D81B60",
    "Mildly integrated" = "#1E88E5"
  )
  
  regime_shapes <- c(
    "Mildly explosive"  = 16,
    "Mildly integrated" = 17
  )
  
  # -----------------------------
  # Camera-ready plot
  # -----------------------------
  
  p <- ggplot(plot_df, aes(x = Item, y = Estimate)) +
    geom_hline(
      aes(yintercept = Reference),
      linewidth = 0.85,
      linetype = "dashed",
      color = "gray25"
    ) +
    geom_errorbar(
      aes(ymin = Lower, ymax = Upper, color = Regime),
      width = 0.18,
      linewidth = 1.25,
      alpha = 0.90
    ) +
    geom_point(
      aes(color = Regime, shape = Regime),
      size = 4.8,
      stroke = 1.2
    ) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
    scale_color_manual(values = regime_cols) +
    scale_shape_manual(values = regime_shapes) +
    labs(
      x = NULL,
      y = NULL,
      color = NULL,
      shape = NULL
    ) +
    theme_bw(base_size = 16) +
    theme(
      strip.background = element_rect(fill = "gray10", color = "gray10"),
      strip.text = element_text(
        face = "bold",
        color = "white",
        size = 17
      ),
      
      axis.text.x = element_text(
        angle = 25,
        hjust = 1,
        face = "bold",
        size = 15,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 15,
        color = "black"
      ),
      
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      axis.ticks.length = grid::unit(0.18, "cm"),
      
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      
      legend.position = "bottom",
      legend.text = element_text(size = 15),
      legend.key.width = grid::unit(1.4, "cm"),
      legend.key.height = grid::unit(0.55, "cm"),
      
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      
      plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 5.2, linewidth = 1.3)
      ),
      shape = guide_legend(
        override.aes = list(size = 5.2, linewidth = 1.3)
      )
    )
  
  print(est_df)
  print(p)
  
  # -----------------------------
  # Save outputs
  # -----------------------------
  
  if (isTRUE(save_tables)) {
    write.csv(est_df, file.path(TABLE_DIR, table_file), row.names = FALSE)
  }
  
  if (isTRUE(save_figure)) {
    ggsave(
      file.path(OUTPUT_DIR, paste0(figure_file_base, ".pdf")),
      plot = p,
      width = 11,
      height = 8
    )
    
    ggsave(
      file.path(OUTPUT_DIR, paste0(figure_file_base, ".png")),
      plot = p,
      width = 11,
      height = 8,
      dpi = 300
    )
  }
  
  invisible(list(
    estimates = est_df,
    plot_data = plot_df,
    plot = p,
    prices = price_df
  ))
}

# ------------------------------------------------------------------------------
# Run: commodity markets, 2022--2026
# ------------------------------------------------------------------------------

out_commodities <- make_delta_gamma_plot(
  asset_df = COMMODITY_ASSETS,
  group_name = "Asset",
  start_date = "2022-01-01",
  end_date = "2026-03-31",
  lambda = 0.05,
  table_file = "commodity_delta_gamma.csv",
  figure_file_base = "commodity_delta_gamma_two_regimes"
)

# ------------------------------------------------------------------------------
# Run: global real-estate markets, 2020--2026
# ------------------------------------------------------------------------------

out_global_re <- make_delta_gamma_plot(
  asset_df = GLOBAL_RE_ASSETS,
  group_name = "Market",
  start_date = "2020-01-01",
  end_date = "2026-03-31",
  lambda = 0.05,
  table_file = "global_real_estate_delta_gamma.csv",
  figure_file_base = "global_real_estate_delta_gamma_two_regimes"
)

