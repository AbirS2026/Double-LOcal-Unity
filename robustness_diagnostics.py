# ============================================================
# Robustness Diagnostics for Moderate-Deviation Simulations
# ============================================================
#
# This script reports two additional robustness checks for the
# finite-sample diagnostics used in the paper.
#
# Robustness Check 1: Non-Gaussian return innovations
# ---------------------------------------------------
# We replace the Gaussian return innovation epsilon_t by a
# standardized Student-t_5 innovation,
#
#     epsilon_t = sqrt(3/5) * Z_t,    Z_t ~ t_5,
#
# so that
#
#     E[epsilon_t] = 0,
#     Var(epsilon_t) = 1,
#     E[epsilon_t^4] < infinity.
#
# The stochastic-volatility innovation eta_t remains Gaussian,
# exactly as in the maintained model.
#
# Robustness Check 2: Alternative initial volatility
# --------------------------------------------------
# We return to Gaussian epsilon_t, but initialize the volatility
# process at sigma_0 = 5 instead of sigma_0 = 1.
#
# In both exercises we report KS diagnostics for:
#
#   (i) the nearly stationary statistic
#
#       T_n = sqrt(n k_n) (rho_hat_n - rho_n),
#
#       with target law N(0, 2c),
#
#   (ii) the mildly explosive statistic
#
#       S_n = rho_n^n k_n (rho_hat_n - rho_n) / (2c),
#
#       with target law standard Cauchy.
#
# For each experiment we compare:
#
#   - Homoskedastic innovations: sigma_t = 1
#   - Stochastic volatility:
#
#       log(sigma_t^2)
#           = phi_n log(sigma_{t-1}^2) + eta_t,
#
#       phi_n = 1 - d / log(log n),
#       eta_t ~ N(0, alpha^2).
#
# Main numerical settings:
#
#       n       = 10,000
#       c       = 0.5
#       d       = 1
#       alpha   = 0.1
#
# The script prints the tables, saves them as CSV files,
# and prints LaTeX versions suitable for the Appendix.
# ============================================================

import numpy as np
import pandas as pd
from scipy.stats import kstest


# ============================================================
# Global settings
# ============================================================

np.random.seed(123)

N = 10_000
C = 0.5
D = 1.0
ALPHA = 0.1

OUTER_REPS = 100
INNER_MC = 100
TOTAL_MC = OUTER_REPS * INNER_MC
BATCH_SIZE = 500

EPS_DF = 5


# ============================================================
# Innovation generators
# ============================================================

def draw_eps_gaussian(size):
    """Gaussian epsilon_t with mean 0 and variance 1."""
    return np.random.randn(size)


def draw_eps_t5(size):
    """
    Standardized Student-t_5 innovation.

    If Z ~ t_nu, then Var(Z) = nu/(nu-2). Therefore,

        epsilon = sqrt((nu-2)/nu) * Z

    has variance 1.

    With nu = 5, epsilon_t is non-Gaussian and has finite
    fourth moment.
    """
    scale = np.sqrt((EPS_DF - 2.0) / EPS_DF)
    return scale * np.random.standard_t(df=EPS_DF, size=size)


# ============================================================
# k_n designs
# ============================================================

K_SPECS_NEAR = {
    r"$10$": lambda n: 10.0,
    r"$\log n$": lambda n: np.log(n),
    r"$n^{0.10}$": lambda n: n**0.10,
    r"$n^{0.25}$": lambda n: n**0.25,
    r"$n^{0.50}$": lambda n: n**0.50,
    r"$n^{0.75}$": lambda n: n**0.75,
    r"$n/\log n$": lambda n: n / np.log(n),
}

K_SPECS_EXPLOSIVE = {
    r"$\log n$": lambda n: np.log(n),
    r"$n^{0.10}$": lambda n: n**0.10,
    r"$n^{0.25}$": lambda n: n**0.25,
    r"$n^{0.50}$": lambda n: n**0.50,
    r"$n^{0.75}$": lambda n: n**0.75,
    r"$n^{0.90}$": lambda n: n**0.90,
    r"$n/\log n$": lambda n: n / np.log(n),
}


# ============================================================
# Core simulation functions
# ============================================================

def simulate_homoskedastic_stats(
    n,
    k,
    c,
    regime,
    eps_draw
):
    """
    Simulate normalized OLS statistics under homoskedasticity.

    Parameters
    ----------
    eps_draw : callable
        Function generating iid epsilon_t innovations.

    Returns
    -------
    ndarray of shape (OUTER_REPS, INNER_MC)
    """

    all_stats = []
    remaining = TOTAL_MC

    while remaining > 0:
        m = min(BATCH_SIZE, remaining)
        remaining -= m

        if regime == "near":
            rho = 1.0 - c / k

            y_prev = np.zeros(m)
            num = np.zeros(m)
            den = np.zeros(m)

            for _ in range(n):
                u = eps_draw(m)

                num += y_prev * u
                den += y_prev**2

                y_prev = rho * y_prev + u

            stats = np.sqrt(n * k) * (num / den)

        elif regime == "explosive":
            rho = 1.0 + c / k
            log_rho = np.log1p(c / k)

            # Scaled state x_t = rho^{-t} y_t
            # avoids direct overflow from rho_n^n.
            x_prev = np.zeros(m)
            num_scaled = np.zeros(m)
            den_scaled = np.zeros(m)

            for t in range(1, n + 1):
                u = eps_draw(m)

                w = np.exp(-(n - t + 1) * log_rho)

                num_scaled += w * x_prev * u
                den_scaled += (w**2) * (x_prev**2)

                x_prev += np.exp(-t * log_rho) * u

            stats = (k / (2.0 * c)) * (num_scaled / den_scaled)

        else:
            raise ValueError("regime must be 'near' or 'explosive'")

        all_stats.append(stats)

    stats = np.concatenate(all_stats)
    return stats.reshape(OUTER_REPS, INNER_MC)


def simulate_sv_stats(
    n,
    k,
    c,
    d,
    alpha,
    regime,
    eps_draw,
    sigma0=1.0
):
    """
    Simulate normalized OLS statistics under stochastic volatility.

    Volatility model:
        z_t = phi_n z_{t-1} + eta_t,
        z_t = log(sigma_t^2),

        phi_n = 1 - d/log(log n),
        eta_t ~ N(0, alpha^2).

    The initial volatility is sigma_0, hence
        z_0 = log(sigma_0^2) = 2 log(sigma_0).

    epsilon_t and eta_t are generated independently.
    """

    if sigma0 <= 0:
        raise ValueError("sigma0 must be strictly positive.")

    phi = 1.0 - d / np.log(np.log(n))

    all_stats = []
    remaining = TOTAL_MC

    while remaining > 0:
        m = min(BATCH_SIZE, remaining)
        remaining -= m

        z = np.full(m, 2.0 * np.log(sigma0))

        if regime == "near":
            rho = 1.0 - c / k

            y_prev = np.zeros(m)
            num = np.zeros(m)
            den = np.zeros(m)

            for _ in range(n):
                eps = eps_draw(m)
                eta = alpha * np.random.randn(m)

                z = phi * z + eta
                u = np.exp(0.5 * z) * eps

                num += y_prev * u
                den += y_prev**2

                y_prev = rho * y_prev + u

            stats = np.sqrt(n * k) * (num / den)

        elif regime == "explosive":
            rho = 1.0 + c / k
            log_rho = np.log1p(c / k)

            x_prev = np.zeros(m)
            num_scaled = np.zeros(m)
            den_scaled = np.zeros(m)

            for t in range(1, n + 1):
                eps = eps_draw(m)
                eta = alpha * np.random.randn(m)

                z = phi * z + eta
                u = np.exp(0.5 * z) * eps

                w = np.exp(-(n - t + 1) * log_rho)

                num_scaled += w * x_prev * u
                den_scaled += (w**2) * (x_prev**2)

                x_prev += np.exp(-t * log_rho) * u

            stats = (k / (2.0 * c)) * (num_scaled / den_scaled)

        else:
            raise ValueError("regime must be 'near' or 'explosive'")

        all_stats.append(stats)

    stats = np.concatenate(all_stats)
    return stats.reshape(OUTER_REPS, INNER_MC)


# ============================================================
# KS diagnostics
# ============================================================

def ks_summary(stats_matrix, regime):
    """
    Run one KS test for each outer replication.

    Returns
    -------
    mean_ks : float
        Mean KS statistic across outer experiments.

    acceptance_prop : float
        Fraction of KS tests with p-value > 0.05.
    """

    ks_vals = []
    accept = []

    for b in range(OUTER_REPS):
        x = stats_matrix[b, :]
        x = x[np.isfinite(x)]

        if len(x) < 2:
            ks_vals.append(np.nan)
            accept.append(np.nan)
            continue

        if regime == "near":
            result = kstest(
                x,
                "norm",
                args=(0.0, np.sqrt(2.0 * C))
            )

        elif regime == "explosive":
            result = kstest(
                x,
                "cauchy",
                args=(0.0, 1.0)
            )

        else:
            raise ValueError("regime must be 'near' or 'explosive'")

        ks_vals.append(result.statistic)
        accept.append(result.pvalue > 0.05)

    return np.nanmean(ks_vals), np.nanmean(accept)


# ============================================================
# Generic runners
# ============================================================

def run_one_design(
    k_name,
    k_value,
    regime,
    model,
    eps_draw,
    sigma0
):
    """
    Run one k_n design under either homoskedasticity or SV.
    """

    if model == "homo":
        stats_matrix = simulate_homoskedastic_stats(
            n=N,
            k=k_value,
            c=C,
            regime=regime,
            eps_draw=eps_draw
        )

    elif model == "sv":
        stats_matrix = simulate_sv_stats(
            n=N,
            k=k_value,
            c=C,
            d=D,
            alpha=ALPHA,
            regime=regime,
            eps_draw=eps_draw,
            sigma0=sigma0
        )

    else:
        raise ValueError("model must be 'homo' or 'sv'")

    mean_ks, acc_prop = ks_summary(stats_matrix, regime)

    return {
        "k_n": k_name,
        "Mean KS stat": mean_ks,
        "Acceptance Proportion": acc_prop
    }


def run_table(
    k_specs,
    regime,
    eps_draw,
    sigma0,
    experiment_name
):
    """
    Run all k_n designs and combine homoskedastic and SV results.
    """

    homo_rows = []
    sv_rows = []

    for k_name, k_fun in k_specs.items():
        k_value = float(k_fun(N))

        print(
            f"{experiment_name} | Homoskedastic | "
            f"regime = {regime:10s} | k_n = {k_name}"
        )

        homo_rows.append(
            run_one_design(
                k_name=k_name,
                k_value=k_value,
                regime=regime,
                model="homo",
                eps_draw=eps_draw,
                sigma0=sigma0
            )
        )

        print(
            f"{experiment_name} | SV             | "
            f"regime = {regime:10s} | k_n = {k_name}"
        )

        sv_rows.append(
            run_one_design(
                k_name=k_name,
                k_value=k_value,
                regime=regime,
                model="sv",
                eps_draw=eps_draw,
                sigma0=sigma0
            )
        )

    homo_df = pd.DataFrame(homo_rows)
    sv_df = pd.DataFrame(sv_rows)

    return homo_df.merge(
        sv_df,
        on="k_n",
        suffixes=("_Homo", "_SV")
    )


# ============================================================
# LaTeX output helper
# ============================================================

def print_latex_combined_table(df, caption, label):
    """
    Print a LaTeX table in the same format as the paper.
    """

    print("\n" + caption)
    print(r"\begin{table}[h!]")
    print(r"\centering")
    print(r"\begin{tabular}{lrrrr}")
    print(r"\toprule")
    print(
        r"& \multicolumn{2}{c}{Homoskedastic} "
        r"& \multicolumn{2}{c}{Stochastic volatility} \\"
    )
    print(r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}")
    print(
        r"\(k_n\) & Mean KS stat & Acceptance prop. "
        r"& Mean KS stat & Acceptance prop. \\"
    )
    print(r"\midrule")

    for _, row in df.iterrows():
        print(
            f"{row['k_n']} & "
            f"{row['Mean KS stat_Homo']:.4f} & "
            f"{row['Acceptance Proportion_Homo']:.2f} & "
            f"{row['Mean KS stat_SV']:.4f} & "
            f"{row['Acceptance Proportion_SV']:.2f} \\\\"
        )

    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(rf"\caption{{{caption}}}")
    print(rf"\label{{{label}}}")
    print(r"\end{table}")


# ============================================================
# ROBUSTNESS CHECK 1:
# Non-Gaussian epsilon_t = standardized Student-t_5
# sigma_0 remains at the baseline value 1
# ============================================================

print("\n")
print("=" * 70)
print("ROBUSTNESS CHECK 1: STANDARDIZED STUDENT-t_5 EPSILON")
print("sigma_0 = 1")
print("=" * 70)

t5_near = run_table(
    K_SPECS_NEAR,
    regime="near",
    eps_draw=draw_eps_t5,
    sigma0=1.0,
    experiment_name="Student-t5"
)

t5_explosive = run_table(
    K_SPECS_EXPLOSIVE,
    regime="explosive",
    eps_draw=draw_eps_t5,
    sigma0=1.0,
    experiment_name="Student-t5"
)

print("\nStudent-t5 robustness: Nearly stationary")
print(t5_near.round(4))

print("\nStudent-t5 robustness: Mildly explosive")
print(t5_explosive.round(4))

t5_near.to_csv(
    "robustness_t5_near_stationary.csv",
    index=False
)

t5_explosive.to_csv(
    "robustness_t5_mild_explosive.csv",
    index=False
)

print_latex_combined_table(
    t5_near,
    caption=(
        r"KS diagnostics for the nearly stationary statistic \(T_n\) "
        r"under standardized Student-\(t_5\) innovations. "
        r"Higher acceptance proportions indicate closer agreement "
        r"with the Gaussian limit \(N(0,2c)\)."
    ),
    label="tab:robust_t5_near"
)

print_latex_combined_table(
    t5_explosive,
    caption=(
        r"KS diagnostics for the mildly explosive statistic \(S_n\) "
        r"under standardized Student-\(t_5\) innovations. "
        r"Higher acceptance proportions indicate closer agreement "
        r"with the standard Cauchy limit."
    ),
    label="tab:robust_t5_explosive"
)


# ============================================================
# ROBUSTNESS CHECK 2:
# Alternative volatility initialization sigma_0 = 5
# epsilon_t returns to Gaussian
# ============================================================

print("\n")
print("=" * 70)
print("ROBUSTNESS CHECK 2: ALTERNATIVE INITIAL VOLATILITY")
print("epsilon_t ~ N(0,1), sigma_0 = 5")
print("=" * 70)

sigma5_near = run_table(
    K_SPECS_NEAR,
    regime="near",
    eps_draw=draw_eps_gaussian,
    sigma0=5.0,
    experiment_name="sigma0=5"
)

sigma5_explosive = run_table(
    K_SPECS_EXPLOSIVE,
    regime="explosive",
    eps_draw=draw_eps_gaussian,
    sigma0=5.0,
    experiment_name="sigma0=5"
)

print("\nAlternative sigma_0 robustness: Nearly stationary")
print(sigma5_near.round(4))

print("\nAlternative sigma_0 robustness: Mildly explosive")
print(sigma5_explosive.round(4))

sigma5_near.to_csv(
    "robustness_sigma0_5_near_stationary.csv",
    index=False
)

sigma5_explosive.to_csv(
    "robustness_sigma0_5_mild_explosive.csv",
    index=False
)

print_latex_combined_table(
    sigma5_near,
    caption=(
        r"KS diagnostics for the nearly stationary statistic \(T_n\) "
        r"under the alternative volatility initialization \(\sigma_0=5\). "
        r"Higher acceptance proportions indicate closer agreement "
        r"with the Gaussian limit \(N(0,2c)\)."
    ),
    label="tab:robust_sigma5_near"
)

print_latex_combined_table(
    sigma5_explosive,
    caption=(
        r"KS diagnostics for the mildly explosive statistic \(S_n\) "
        r"under the alternative volatility initialization \(\sigma_0=5\). "
        r"Higher acceptance proportions indicate closer agreement "
        r"with the standard Cauchy limit."
    ),
    label="tab:robust_sigma5_explosive"
)


# ============================================================
# End
# ============================================================

print("\nAll robustness diagnostics completed.")
