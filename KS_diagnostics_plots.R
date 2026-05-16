# ============================================================
#  diagnostic plot for Figure 1
# Gaussian and Cauchy limits under stochastic volatility
# ============================================================

library(ggplot2)
library(patchwork)
library(cowplot)

set.seed(1233)

# -----------------------------
# Helper functions
# -----------------------------

simulate_sv_innovations <- function(n, d, alpha) {
  phi <- 1 - d / log(log(n))
  z <- numeric(n + 1)
  eps <- rnorm(n)
  eta <- rnorm(n, mean = 0, sd = alpha)
  
  u <- numeric(n)
  
  for (t in 1:n) {
    z[t + 1] <- phi * z[t] + eta[t]
    sigma_t <- exp(0.5 * z[t + 1])
    u[t] <- sigma_t * eps[t]
  }
  
  u
}

simulate_near_stat <- function(n, c, d, alpha, k) {
  rho <- 1 - c / k
  u <- simulate_sv_innovations(n, d, alpha)
  
  y <- numeric(n + 1)
  
  for (t in 1:n) {
    y[t + 1] <- rho * y[t] + u[t]
  }
  
  y_lag <- y[1:n]
  rho_hat <- rho + sum(y_lag * u) / sum(y_lag^2)
  
  sqrt(n * k) * (rho_hat - rho)
}

simulate_explosive_stat <- function(n, c, d, alpha, k) {
  rho <- 1 + c / k
  log_rho <- log(rho)
  u <- simulate_sv_innovations(n, d, alpha)
  
  # Scaled recursion: x_t = rho^{-t} y_t
  x_prev <- 0
  num_scaled <- 0
  den_scaled <- 0
  
  for (t in 1:n) {
    w <- exp(-(n - t + 1) * log_rho)
    
    num_scaled <- num_scaled + w * x_prev * u[t]
    den_scaled <- den_scaled + w^2 * x_prev^2
    
    x_prev <- x_prev + exp(-t * log_rho) * u[t]
  }
  
  (k / (2 * c)) * num_scaled / den_scaled
}

# -----------------------------
# Specifications
# -----------------------------

B <- 5000

beta_near <- 0.1
beta_exp  <- 0.1

# Nearly stationary
n_near <- 1000
c_near <- 0.5
d_near <- 1
alpha_near <- 0.5
k_near <- n_near^beta_near

# Mildly explosive
n_exp <- 1000
c_exp <- 0.5
d_exp <- 1
alpha_exp <- 0.5
k_exp <- n_exp^beta_exp

# -----------------------------
# Simulate statistics
# -----------------------------

T_vals <- replicate(
  B,
  simulate_near_stat(
    n = n_near,
    c = c_near,
    d = d_near,
    alpha = alpha_near,
    k = k_near
  )
)

S_vals <- replicate(
  B,
  simulate_explosive_stat(
    n = n_exp,
    c = c_exp,
    d = d_exp,
    alpha = alpha_exp,
    k = k_exp
  )
)

# KS diagnostics, using full samples
ks_T <- ks.test(T_vals, "pnorm", mean = 0, sd = sqrt(2 * c_near))
ks_S <- ks.test(S_vals, "pcauchy", location = 0, scale = 1)

cat("Nearly stationary KS:",
    round(as.numeric(ks_T$statistic), 4),
    "p-value:",
    round(ks_T$p.value, 4), "\n")

cat("Mildly explosive KS:",
    round(as.numeric(ks_S$statistic), 4),
    "p-value:",
    round(ks_S$p.value, 4), "\n")

df_T <- data.frame(statistic = T_vals)
df_S <- data.frame(statistic = S_vals)

# Plotting-only clipping for Cauchy panel.
# The KS test above still uses the full S_vals sample.
x_low_S <- quantile(S_vals, 0.005, na.rm = TRUE)
x_high_S <- quantile(S_vals, 0.995, na.rm = TRUE)
df_S_plot <- subset(df_S, statistic >= x_low_S & statistic <= x_high_S)

# -----------------------------
# Camera-ready theme
# -----------------------------

base_theme <- theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 17, color = "black"),
    axis.text.y = element_text(size = 17, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = "grey80"),
    plot.margin = margin(t = 10, r = 14, b = 10, l = 14)
  )

# -----------------------------
# Left panel: Gaussian limit
# -----------------------------

p_normal <- ggplot(df_T, aes(x = statistic)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 35,
    fill = "steelblue",
    alpha = 0.45,
    color = "white"
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = 0, sd = sqrt(2 * c_near)),
    color = "red",
    linewidth = 1.25
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(x = NULL, y = NULL) +
  base_theme

# -----------------------------
# Right panel: Cauchy limit
# -----------------------------

p_cauchy <- ggplot(df_S_plot, aes(x = statistic)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 60,
    fill = "steelblue",
    alpha = 0.35,
    color = "white"
  ) +
  stat_function(
    fun = dcauchy,
    args = list(location = 0, scale = 1),
    color = "red",
    linewidth = 1.25,
    n = 1000,
    xlim = c(x_low_S, x_high_S)
  ) +
  coord_cartesian(xlim = c(x_low_S, x_high_S)) +
  labs(x = NULL, y = NULL) +
  base_theme

# -----------------------------
# Combine panels with outside y-label
# -----------------------------

combined_plot <- p_normal + p_cauchy + plot_layout(nrow = 1)

final_plot <- ggdraw(combined_plot) +
  draw_label(
    "Density",
    x = 0.01,      # farther from the panels
    y = 0.5,
    angle = 90,
    size = 18,
    fontface = "bold"
  )

print(final_plot)

# -----------------------------
# Save figure
# -----------------------------

ggsave(
  filename = paste0(
    "diagnostic_normal_cauchy_beta_near_", beta_near,
    "_beta_exp_", beta_exp, ".pdf"
  ),
  plot = final_plot,
  width = 11.5,
  height = 4.6
)

ggsave(
  filename = paste0(
    "diagnostic_normal_cauchy_beta_near_", beta_near,
    "_beta_exp_", beta_exp, ".png"
  ),
  plot = final_plot,
  width = 11.5,
  height = 4.6,
  dpi = 300
)