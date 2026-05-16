# ============================================================
# Camera-ready volatility motivation figure
# Figure 2: price paths and rolling volatility
#
# Output:
#   figures/volatility_motivation_2x2_price_rolling_sd.pdf
#   figures/volatility_motivation_2x2_price_rolling_sd.png
# ============================================================

library(quantmod)
library(dplyr)
library(ggplot2)
library(zoo)
library(gridExtra)
library(grid)

# -----------------------------
# Helper functions
# -----------------------------

fetch_price_series <- function(ticker, start_date, end_date) {
  x <- quantmod::getSymbols(
    ticker,
    src = "yahoo",
    from = as.Date(start_date),
    to = as.Date(end_date),
    auto.assign = FALSE,
    warnings = FALSE
  )
  
  px <- tryCatch(quantmod::Ad(x), error = function(e) quantmod::Cl(x))
  na.omit(px)
}

compute_rolling_sd <- function(px_xts, k_lags = 50) {
  zoo::rollapply(
    px_xts,
    width = k_lags,
    FUN = sd,
    align = "right",
    fill = NA
  )
}

theme_paper <- function(base_size = 15) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_blank()
    )
}

# -----------------------------
# Figure specifications
# -----------------------------

PAPER_END_DATE <- as.Date("2026-04-30")

fig2_assets <- data.frame(
  Asset = c(
    "NVIDIA: AI cycle",
    "Natural gas: subprime-crisis cycle",
    "Bitcoin: crypto (2020-2021) cycle",
    "Cocoa: recent commodity cycle"
  ),
  Ticker = c("NVDA", "NG=F", "BTC-USD", "CC=F"),
  Start_Date = as.Date(c(
    "2022-01-01",
    "2006-01-01",
    "2020-01-01",
    "2023-09-01"
  )),
  End_Date = as.Date(c(
    as.character(PAPER_END_DATE),
    "2009-12-31",
    "2022-12-31",
    as.character(PAPER_END_DATE)
  )),
  stringsAsFactors = FALSE
)

fig2_k_lags <- 50
target_fraction <- 0.50

stock_color <- "#1f77b4"
volatility_color <- "#d62728"

# -----------------------------
# Panel function
# -----------------------------

make_price_vol_panel <- function(asset_name, ticker, start_date, end_date) {
  
  px_xts <- fetch_price_series(ticker, start_date, end_date)
  
  if (is.null(px_xts) || NROW(px_xts) == 0) {
    stop(paste("Failed to load", ticker, "from Yahoo Finance."))
  }
  
  sd_xts <- compute_rolling_sd(px_xts, k_lags = fig2_k_lags)
  colnames(sd_xts) <- "RollingSD"
  
  df <- data.frame(
    Date = as.Date(index(px_xts)),
    Price = as.numeric(px_xts),
    stringsAsFactors = FALSE
  ) |>
    left_join(
      data.frame(
        Date = as.Date(index(sd_xts)),
        RollingSD = as.numeric(sd_xts),
        stringsAsFactors = FALSE
      ),
      by = "Date"
    )
  
  price_top <- quantile(df$Price, 0.98, na.rm = TRUE)
  sd_top <- quantile(df$RollingSD, 0.98, na.rm = TRUE)
  
  scale_factor <- if (is.finite(sd_top) && sd_top > 0) {
    target_fraction * price_top / sd_top
  } else {
    1
  }
  
  df$RollingSD_scaled <- df$RollingSD * scale_factor
  
  ggplot(df, aes(x = Date)) +
    geom_line(
      aes(y = Price),
      color = stock_color,
      linewidth = 1.05
    ) +
    geom_line(
      aes(y = RollingSD_scaled),
      color = volatility_color,
      linewidth = 1.05,
      alpha = 0.90,
      na.rm = TRUE
    ) +
    scale_y_continuous(
      name = NULL,
      sec.axis = sec_axis(~ . / scale_factor, name = NULL)
    ) +
    scale_x_date(
      date_breaks = "6 months",
      date_labels = "%Y %b"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = asset_name
    ) +
    theme_paper(base_size = 15) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      
      axis.text.y.left = element_text(
        size = 15,
        color = stock_color,
        face = "bold"
      ),
      axis.text.y.right = element_text(
        size = 15,
        color = volatility_color,
        face = "bold"
      ),
      axis.ticks.y.left = element_line(color = stock_color, linewidth = 0.7),
      axis.ticks.y.right = element_line(color = volatility_color, linewidth = 0.7),
      
      axis.text.x = element_text(
        size = 13,
        angle = 45,
        hjust = 1,
        color = "black"
      ),
      
      plot.margin = margin(t = 8, r = 12, b = 8, l = 12),
      panel.grid.minor = element_blank()
    )
}

# -----------------------------
# Build panels
# -----------------------------

fig2_panels <- lapply(seq_len(nrow(fig2_assets)), function(i) {
  make_price_vol_panel(
    asset_name = fig2_assets$Asset[i],
    ticker = fig2_assets$Ticker[i],
    start_date = fig2_assets$Start_Date[i],
    end_date = fig2_assets$End_Date[i]
  )
})

# -----------------------------
# Align all four panels exactly
# This fixes unequal box sizes caused by different axis-label widths/heights.
# -----------------------------

fig2_grobs <- lapply(fig2_panels, ggplotGrob)

max_widths <- do.call(grid::unit.pmax, lapply(fig2_grobs, function(g) g$widths))
max_heights <- do.call(grid::unit.pmax, lapply(fig2_grobs, function(g) g$heights))

fig2_grobs <- lapply(fig2_grobs, function(g) {
  g$widths <- max_widths
  g$heights <- max_heights
  g
})

middle_panel <- gridExtra::arrangeGrob(
  grobs = fig2_grobs,
  nrow = 2,
  ncol = 2,
  widths = c(1, 1),
  heights = c(1, 1)
)

# -----------------------------
# Outer axis labels
# -----------------------------

left_axis_label <- grid::textGrob(
  "Asset Price",
  rot = 90,
  gp = grid::gpar(
    col = stock_color,
    fontsize = 22,
    fontface = "bold"
  )
)

right_axis_label <- grid::textGrob(
  "Price Volatility",
  rot = 270,
  gp = grid::gpar(
    col = volatility_color,
    fontsize = 22,
    fontface = "bold"
  )
)

fig2 <- gridExtra::arrangeGrob(
  grobs = list(left_axis_label, middle_panel, right_axis_label),
  ncol = 3,
  widths = c(0.07, 0.86, 0.07)
)

# -----------------------------
# Save figure
# -----------------------------

if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

grDevices::pdf(
  file = "figures/volatility_motivation_2x2_price_rolling_sd.pdf",
  width = 13.5,
  height = 8.8
)
grid::grid.newpage()
grid::grid.draw(fig2)
grDevices::dev.off()

grDevices::png(
  filename = "figures/volatility_motivation_2x2_price_rolling_sd.png",
  width = 13.5,
  height = 8.8,
  units = "in",
  res = 300
)
grid::grid.newpage()
grid::grid.draw(fig2)
grDevices::dev.off()

# Display in RStudio / notebook
grid::grid.newpage()
grid::grid.draw(fig2)

