# ============================================================
# Replication code for main-paper Figure 4 and Table 3
#
# Figure 4:
#   Historical bubble candidates with detected exuberance windows.
#   Output:
#     figures/historical_bubble_candidates_2x2.pdf
#     figures/historical_bubble_candidates_2x2.png
#
# Table 3:
#   Comparison of the stochastic-volatility robust diagnostic
#   with the homoskedastic PWY-type diagnostic across canonical
#   bubble windows and recent comparison windows.
#   Output:
#     tables/sv_adf_vs_pwy_canonical_vs_recent.csv
#     tables/table3_canonical_recent_comparison.csv
#     tables/table3_canonical_recent_comparison.tex
# ============================================================

library(quantmod)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# ------------------------------------------------------------
# Global settings
# ------------------------------------------------------------

SAVE_FIGURES <- TRUE
SAVE_TABLES  <- TRUE

OUTPUT_DIR <- "figures"
TABLE_DIR  <- "tables"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(TABLE_DIR)) dir.create(TABLE_DIR, recursive = TRUE)

# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------

FIG4_PDF <- file.path(OUTPUT_DIR, "historical_bubble_candidates_2x2.pdf")
FIG4_PNG <- file.path(OUTPUT_DIR, "historical_bubble_candidates_2x2.png")

FIG4_TABLE_CSV <- file.path(TABLE_DIR, "historical_bubble_candidates_results.csv")

TABLE3_FULL_CSV <- file.path(TABLE_DIR, "sv_adf_vs_pwy_canonical_vs_recent.csv")
TABLE3_CSV      <- file.path(TABLE_DIR, "table3_canonical_recent_comparison.csv")
TABLE3_TEX      <- file.path(TABLE_DIR, "table3_canonical_recent_comparison.tex")

# Optional: if FRED access fails, put CSUSHPINSA.csv in your working directory.
HPI_CSV <- "CSUSHPINSA.csv"

# ============================================================
# Part I: Figure 4 historical bubble candidates
# ============================================================

# ------------------------------------------------------------
# Episode definitions for Figure 4
# ------------------------------------------------------------

episode_assets <- data.frame(
  Asset = c(
    "U.S. House Price Index",
    "Crude Oil",
    "Nasdaq",
    "Bitcoin"
  ),
  Ticker = c(
    "CSUSHPINSA",
    "CL=F",
    "^IXIC",
    "BTC-USD"
  ),
  Source = c(
    "FRED_HPI",
    "YAHOO",
    "YAHOO",
    "YAHOO"
  ),
  Plot_Start = as.Date(c(
    "2000-01-01",
    "2000-01-01",
    "1995-01-01",
    "2020-01-01"
  )),
  Plot_End = as.Date(c(
    "2008-12-31",
    "2008-12-31",
    "2002-12-31",
    "2022-12-31"
  )),
  Test_Start = as.Date(c(
    "2002-01-01",
    "2007-01-01",
    "1995-01-01",
    "2020-03-15"
  )),
  Test_End = as.Date(c(
    "2006-12-31",
    "2008-08-01",
    "2000-09-30",
    "2021-11-10"
  )),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Data helpers for Figure 4
# ------------------------------------------------------------

get_yahoo_px <- function(ticker, asset, start, end, min_obs = 40) {
  x <- tryCatch(
    quantmod::getSymbols(
      ticker,
      src = "yahoo",
      from = as.Date(start),
      to = as.Date(end),
      auto.assign = FALSE,
      warnings = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(x)) {
    warning(paste("Could not fetch:", asset, ticker))
    return(NULL)
  }
  
  px <- tryCatch(quantmod::Ad(x), error = function(e) quantmod::Cl(x))
  
  out <- data.frame(
    Date = as.Date(index(px)),
    Asset = asset,
    Price = as.numeric(px),
    stringsAsFactors = FALSE
  ) |>
    na.omit()
  
  if (nrow(out) < min_obs) {
    warning(paste("Too few observations:", asset, ticker))
    return(NULL)
  }
  
  out
}

get_hpi_px <- function(asset, start, end, local_file = HPI_CSV) {
  x <- tryCatch(
    quantmod::getSymbols(
      "CSUSHPINSA",
      src = "FRED",
      from = as.Date(start),
      to = as.Date(end),
      auto.assign = FALSE,
      warnings = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(x)) {
    return(
      data.frame(
        Date = as.Date(index(x)),
        Asset = asset,
        Price = as.numeric(x[, 1]),
        stringsAsFactors = FALSE
      ) |>
        na.omit()
    )
  }
  
  if (!file.exists(local_file)) {
    stop(
      paste0(
        "FRED failed and local file not found: ", local_file,
        ". Download CSUSHPINSA.csv and place it in getwd()."
      )
    )
  }
  
  z <- read.csv(local_file)
  names(z)[1:2] <- c("Date", "Price")
  
  z |>
    mutate(
      Date = as.Date(Date),
      Price = as.numeric(Price),
      Asset = asset
    ) |>
    filter(Date >= as.Date(start), Date <= as.Date(end)) |>
    select(Date, Asset, Price) |>
    na.omit()
}

fetch_episode_px <- function(asset, ticker, source, start, end) {
  if (source == "FRED_HPI") {
    get_hpi_px(asset, start, end)
  } else {
    get_yahoo_px(ticker, asset, start, end)
  }
}

# ------------------------------------------------------------
# Basic DF-type statistic used in diagnostics
# ------------------------------------------------------------

compute_df_delta <- function(X, normalize_start = FALSE) {
  X <- as.numeric(X)
  X <- X[is.finite(X)]
  
  if (length(X) < 3) {
    return(list(n_eff = NA_real_, delta_hat = NA_real_, df_delta = NA_real_))
  }
  
  if (normalize_start) X <- X / X[1]
  
  y <- X[-1]
  xlag <- X[-length(X)]
  n_eff <- length(y)
  
  denom <- sum(xlag^2)
  if (!is.finite(denom) || denom <= 0) {
    return(list(n_eff = n_eff, delta_hat = NA_real_, df_delta = NA_real_))
  }
  
  delta_hat <- sum(y * xlag) / denom
  df_delta <- n_eff * (delta_hat - 1)
  
  list(
    n_eff = n_eff,
    delta_hat = delta_hat,
    df_delta = df_delta
  )
}

# ------------------------------------------------------------
# Create Figure 4
# ------------------------------------------------------------

make_historical_bubble_plot <- function(
    assets = episode_assets,
    normalize_start = TRUE
) {
  
  price_df <- bind_rows(lapply(seq_len(nrow(assets)), function(i) {
    fetch_episode_px(
      asset = assets$Asset[i],
      ticker = assets$Ticker[i],
      source = assets$Source[i],
      start = assets$Plot_Start[i],
      end = assets$Plot_End[i]
    )
  }))
  
  # Compute diagnostic over the marked red-window period.
  results <- bind_rows(lapply(seq_len(nrow(assets)), function(i) {
    
    asset_i <- assets$Asset[i]
    start_i <- assets$Test_Start[i]
    end_i   <- assets$Test_End[i]
    
    tmp <- price_df |>
      filter(Asset == asset_i, Date >= start_i, Date <= end_i) |>
      arrange(Date)
    
    X <- tmp$Price
    if (normalize_start && length(X) > 0) X <- X / X[1]
    
    stat <- compute_df_delta(X)
    
    data.frame(
      Asset = asset_i,
      Test_Start = as.character(start_i),
      Test_End = as.character(end_i),
      N_Obs = nrow(tmp),
      N_Eff = stat$n_eff,
      Delta_Hat = stat$delta_hat,
      DF_Delta = stat$df_delta,
      Cutoff_SV = log(stat$n_eff) / 10,
      Cutoff_PWY = -0.08,
      Decision_SV = ifelse(
        stat$df_delta > log(stat$n_eff) / 10,
        "Bubble",
        "No bubble"
      ),
      Decision_PWY = ifelse(
        stat$df_delta > -0.08,
        "Bubble",
        "No bubble"
      ),
      stringsAsFactors = FALSE
    )
  }))
  
  # Normalize each displayed price path to start at 100.
  plot_df <- price_df |>
    left_join(
      assets |>
        select(Asset, Plot_Start, Plot_End, Test_Start, Test_End),
      by = "Asset"
    ) |>
    group_by(Asset) |>
    arrange(Date, .by_group = TRUE) |>
    mutate(Index = 100 * Price / first(Price)) |>
    ungroup()
  
  plot_df$Asset <- factor(plot_df$Asset, levels = assets$Asset)
  
  p <- ggplot(plot_df, aes(x = Date, y = Index)) +
    geom_line(
      linewidth = 1.05,
      color = "#0047AB"
    ) +
    geom_vline(
      aes(xintercept = Test_Start),
      color = "#d62728",
      linewidth = 0.95,
      linetype = "solid"
    ) +
    geom_vline(
      aes(xintercept = Test_End),
      color = "#d62728",
      linewidth = 0.95,
      linetype = "solid"
    ) +
    facet_wrap(~ Asset, scales = "free", nrow = 2, ncol = 2) +
    labs(
      x = NULL,
      y = "Price Index"
    ) +
    theme_bw(base_size = 16) +
    theme(
      strip.background = element_rect(fill = "gray10", color = "gray10"),
      strip.text = element_text(
        face = "bold",
        color = "white",
        size = 18
      ),
      axis.text.x = element_text(
        angle = 35,
        hjust = 1,
        size = 15,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 15,
        color = "black"
      ),
      axis.title.y = element_text(
        size = 19,
        face = "bold",
        margin = margin(r = 14)
      ),
      axis.ticks = element_line(linewidth = 0.65, color = "black"),
      axis.ticks.length = grid::unit(0.18, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
      plot.caption = element_blank(),
      plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
    )
  
  print(results)
  print(p)
  
  if (isTRUE(SAVE_TABLES)) {
    write.csv(results, FIG4_TABLE_CSV, row.names = FALSE)
  }
  
  if (isTRUE(SAVE_FIGURES)) {
    ggsave(FIG4_PDF, plot = p, width = 12.5, height = 8.8)
    ggsave(FIG4_PNG, plot = p, width = 12.5, height = 8.8, dpi = 300)
  }
  
  invisible(list(
    results = results,
    plot = p,
    prices = price_df,
    plot_data = plot_df
  ))
}

# ============================================================
# Part II: Table 3 canonical-vs-recent comparison
# ============================================================

# ------------------------------------------------------------
# Canonical bubble windows and recent comparison windows
# ------------------------------------------------------------

ASSET_PERIODS <- data.frame(
  Asset = c(
    "Nikkei 225",
    "Nasdaq Composite",
    "Shanghai Composite",
    "Crude Oil",
    "Real Estate ETF",
    "Bitcoin",
    "Ethereum",
    "GameStop",
    "ARK Innovation ETF"
  ),
  Ticker = c(
    "^N225",
    "^IXIC",
    "000001.SS",
    "CL=F",
    "IYR",
    "BTC-USD",
    "ETH-USD",
    "GME",
    "ARKK"
  ),
  Bubble_Start = as.Date(c(
    "1985-01-01",
    "1995-01-01",
    "2005-01-01",
    "2007-01-01",
    "2004-01-01",
    "2020-03-15",
    "2020-03-15",
    "2020-12-01",
    "2020-03-15"
  )),
  Bubble_End = as.Date(c(
    "1989-12-29",
    "2000-03-10",
    "2007-10-16",
    "2008-07-11",
    "2007-02-01",
    "2021-11-10",
    "2021-11-10",
    "2021-02-25",
    "2021-02-12"
  )),
  Recent_Start = as.Date(c(
    "2024-01-01",
    "2020-01-01",
    "2020-01-01",
    "2020-01-01",
    "2020-01-01",
    "2024-01-01",
    "2024-01-01",
    "2025-01-01",
    "2025-01-01"
  )),
  Recent_End = as.Date(c(
    "2026-03-31",
    "2026-03-31",
    "2026-03-31",
    "2026-03-31",
    "2026-03-31",
    "2026-03-31",
    "2026-03-31",
    "2026-03-31",
    "2026-03-31"
  )),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Run one diagnostic window
# ------------------------------------------------------------

run_one_test <- function(asset, ticker, period_label, start_date, end_date) {
  px <- get_yahoo_px(
    ticker = ticker,
    asset = asset,
    start = start_date,
    end = end_date,
    min_obs = 40
  )
  
  if (is.null(px)) {
    return(data.frame(
      Asset = asset,
      Ticker = ticker,
      Period = period_label,
      Start_Date = as.character(start_date),
      End_Date = as.character(end_date),
      N_Obs = NA_integer_,
      N_Eff = NA_real_,
      Delta_Hat = NA_real_,
      DF_Delta = NA_real_,
      Cutoff_SV = NA_real_,
      Cutoff_PWY = -0.08,
      Margin_SV = NA_real_,
      Margin_PWY = NA_real_,
      Decision_SV = NA_character_,
      Decision_PWY = NA_character_,
      Agreement = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  
  stat <- compute_df_delta(px$Price, normalize_start = TRUE)
  
  cutoff_sv <- log(stat$n_eff) / 10
  cutoff_pwy <- -0.08
  
  decision_sv <- ifelse(
    is.finite(stat$df_delta) && stat$df_delta > cutoff_sv,
    "Bubble",
    "No bubble"
  )
  
  decision_pwy <- ifelse(
    is.finite(stat$df_delta) && stat$df_delta > cutoff_pwy,
    "Bubble",
    "No bubble"
  )
  
  data.frame(
    Asset = asset,
    Ticker = ticker,
    Period = period_label,
    Start_Date = as.character(start_date),
    End_Date = as.character(end_date),
    N_Obs = nrow(px),
    N_Eff = stat$n_eff,
    Delta_Hat = stat$delta_hat,
    DF_Delta = stat$df_delta,
    Cutoff_SV = cutoff_sv,
    Cutoff_PWY = cutoff_pwy,
    Margin_SV = stat$df_delta - cutoff_sv,
    Margin_PWY = stat$df_delta - cutoff_pwy,
    Decision_SV = decision_sv,
    Decision_PWY = decision_pwy,
    Agreement = ifelse(decision_sv == decision_pwy, "Agree", "Disagree"),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Build full comparison and Table 3
# ------------------------------------------------------------

run_canonical_vs_recent <- function(asset_periods = ASSET_PERIODS) {
  
  # Full diagnostic output: both periods for each asset.
  results <- bind_rows(lapply(seq_len(nrow(asset_periods)), function(i) {
    
    a <- asset_periods$Asset[i]
    t <- asset_periods$Ticker[i]
    
    bind_rows(
      run_one_test(
        asset = a,
        ticker = t,
        period_label = "Canonical episode",
        start_date = asset_periods$Bubble_Start[i],
        end_date = asset_periods$Bubble_End[i]
      ),
      run_one_test(
        asset = a,
        ticker = t,
        period_label = "Recent comparison",
        start_date = asset_periods$Recent_Start[i],
        end_date = asset_periods$Recent_End[i]
      )
    )
  }))
  
  results$Asset <- factor(results$Asset, levels = asset_periods$Asset)
  results$Period <- factor(
    results$Period,
    levels = c("Canonical episode", "Recent comparison")
  )
  
  # Main-paper Table 3: compact classification table.
  table3 <- results |>
    select(Asset, Period, Decision_SV, Decision_PWY) |>
    pivot_wider(
      names_from = Period,
      values_from = c(Decision_SV, Decision_PWY)
    ) |>
    transmute(
      Asset,
      Canonical_SV = `Decision_SV_Canonical episode`,
      Canonical_PWY = `Decision_PWY_Canonical episode`,
      Recent_SV = `Decision_SV_Recent comparison`,
      Recent_PWY = `Decision_PWY_Recent comparison`
    )
  
  # Print full diagnostics and compact Table 3.
  print(results)
  print(table3)
  
  if (isTRUE(SAVE_TABLES)) {
    write.csv(results, TABLE3_FULL_CSV, row.names = FALSE)
    write.csv(table3, TABLE3_CSV, row.names = FALSE)
    
    # LaTeX version of Table 3.
    con <- file(TABLE3_TEX, open = "wt")
    
    writeLines("\\begin{table}[H]", con)
    writeLines("\\centering", con)
    writeLines("\\begin{tabular}{lcccc}", con)
    writeLines("\\toprule", con)
    writeLines("& \\multicolumn{2}{c}{Canonical episode} & \\multicolumn{2}{c}{Recent comparison} \\\\", con)
    writeLines("\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}", con)
    writeLines("Asset & SV-robust & PWY & SV-robust & PWY \\\\", con)
    writeLines("\\midrule", con)
    
    for (i in seq_len(nrow(table3))) {
      writeLines(
        paste0(
          table3$Asset[i], " & ",
          table3$Canonical_SV[i], " & ",
          table3$Canonical_PWY[i], " & ",
          table3$Recent_SV[i], " & ",
          table3$Recent_PWY[i], " \\\\"
        ),
        con
      )
    }
    
    writeLines("\\bottomrule", con)
    writeLines("\\end{tabular}", con)
    writeLines("\\caption{Comparison of bubble classifications from the stochastic-volatility robust diagnostic and the homoskedastic PWY-type diagnostic. The canonical windows correspond to historically recognized exuberance episodes, while the recent comparison windows apply the same diagnostics to recent data for the same assets.}", con)
    writeLines("\\label{tab:canonical_recent_comparison}", con)
    writeLines("\\end{table}", con)
    
    close(con)
  }
  
  invisible(list(
    results = results,
    table3 = table3
  ))
}

# ============================================================
# Run both main-paper outputs
# ============================================================

# Figure 4
out_historical_bubbles <- make_historical_bubble_plot(
  assets = episode_assets,
  normalize_start = TRUE
)

# Table 3
out_compare <- run_canonical_vs_recent()

options(pillar.sigfig = 8)
