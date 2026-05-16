# Double Local-to-Unity: Estimation under Nearly Nonstationary Volatility

This repository contains the replication code for the paper:

**Double Local-to-Unity: Estimation under Nearly Nonstationary Volatility**  
Abir Sarkar and Martin T. Wells  
Paper: [arXiv:2512.06823](https://arxiv.org/pdf/2512.06823)

The paper develops volatility-robust moderate-deviation inference for autoregressive roots when both the conditional mean and stochastic-volatility dynamics are highly persistent. The code in this repository reproduces the simulation diagnostics, volatility-motivation figures, empirical bubble-growth diagnostics, and comparisons with homoskedastic PWY-type bubble-detection procedures.

---

## Repository Structure

### 1. Monte Carlo KS diagnostics

File: [`Tables_1_2_KS_diagnostics.ipynb`](https://github.com/AbirS2026/Double-LOcal-Unity/blob/main/Tables_1_2_KS_diagnostics.ipynb)

This notebook generates the Monte Carlo Kolmogorov--Smirnov diagnostics reported in **Tables 1 and 2** of the paper.

It simulates AR(1) processes under both homoskedastic innovations and nearly nonstationary stochastic-volatility innovations. The notebook compares the simulated normalized OLS statistics with their theoretical limits: the Gaussian limit in the nearly stationary regime and the Cauchy limit in the mildly explosive regime.

Relevant paper outputs:

- Table 1: KS diagnostics for the nearly stationary statistic
- Table 2: KS diagnostics for the mildly explosive statistic



---

### 2. Monte Carlo distributional diagnostic plots

File: [`KS_diagnostics_plots.R`](https://github.com/AbirS2026/Double-LOcal-Unity/blob/main/KS_diagnostics_plots.R)

This script generates the visual Monte Carlo diagnostic plot reported as **Figure 1** of the paper.

The script simulates the normalized OLS statistics under stochastic volatility and overlays their empirical histograms with the corresponding limiting densities: \(N(0,2c)\) for the nearly stationary statistic and the standard Cauchy density for the mildly explosive statistic.

Relevant paper output:

- Figure 1



---

### 3. Volatility motivation across exuberance cycles

File: [`Volatility_spikes_cycle.R`](https://github.com/AbirS2026/Double-LOcal-Unity/blob/main/Volatility_spikes_cycle.R)

This script generates the volatility-motivation figure reported as **Figure 2** of the paper.

The script downloads price data and computes rolling volatility for four episodes of market exuberance or stress:

- NVIDIA during the recent AI cycle
- Natural gas around the 2007--09 crisis period
- Bitcoin during the 2020--21 crypto run-up
- Cocoa during the recent commodity surge

The purpose is to show that volatility is not well described by a fixed-scale homoskedastic specification during periods of exuberance and stress.

Relevant paper output:

- Figure 2



---

### 4. Bubble-growth diagnostics for recent global markets

File: [`Different_global_markets_bubble_growth_rate.R`](https://github.com/AbirS2026/Double-LOcal-Unity/blob/main/Different_global_markets_bubble_growth_rate.R)

This script generates the empirical \((\widehat\rho_n,\widehat\gamma_n)\)-based diagnostics reported in **Figure 3** of the paper.

The script applies the proposed bubble-detection diagnostics to commodity markets over 2022--2026 and global real-estate markets over 2020--2026. It reports estimated autoregressive coefficients, local-deviation parameters, asymptotic confidence intervals, and regime classifications.

Relevant paper output:

- Figure 3



---

### 5. Historical bubble candidates and PWY comparison

File: [`Historical_Comparison_PWY.R`](https://github.com/AbirS2026/Double-LOcal-Unity/blob/main/Historical_Comparison_PWY.R)

This script generates the historical bubble diagnostics and the comparison between the proposed stochastic-volatility robust diagnostic and a homoskedastic PWY-type rule.

The first part of the script produces **Figure 4**, which plots historical bubble candidates and the detected price-exuberance windows. The historical examples include the U.S. housing boom, the crude-oil run-up around the subprime crisis, the late-1990s Nasdaq technology bubble, and Bitcoin during the 2020--21 crypto cycle.

The second part of the script produces **Table 3**, which compares bubble classifications from the stochastic-volatility robust diagnostic and the homoskedastic PWY-type diagnostic across canonical bubble episodes and recent comparison windows.

Relevant paper outputs:

- Figure 4
- Table 3



---

## How to Use the Code

Each script can be run independently. A typical workflow is:

```r
source("KS_diagnostics_plots.R")
source("Volatility_spikes_cycle.R")
source("Different_global_markets_bubble_growth_rate.R")
source("Historical_Comparison_PWY.R")
```

The Jupyter notebook can be run directly:

```text
Tables_1_2_KS_diagnostics.ipynb
```

Suggested order for reproducing the paper outputs:

```text
1. Tables_1_2_KS_diagnostics.ipynb
2. KS_diagnostics_plots.R
3. Volatility_spikes_cycle.R
4. Different_global_markets_bubble_growth_rate.R
5. Historical_Comparison_PWY.R
```

The scripts automatically create the following folders when needed:

```text
figures/
tables/
```

---

## Required Packages

### R packages

```r
install.packages(c(
  "quantmod",
  "dplyr",
  "tidyr",
  "ggplot2",
  "zoo",
  "gridExtra",
  "patchwork",
  "cowplot"
))
```

The `grid` package is also used in some plotting routines, but it is included with base R and usually does not need to be installed separately.

### Python packages

```python
numpy
pandas
scipy
```

---

## Data Sources

Most empirical scripts download daily adjusted closing prices directly from Yahoo Finance using the `quantmod` package:

[https://finance.yahoo.com](https://finance.yahoo.com)

The U.S. house price index is obtained from FRED when available. If FRED access fails, place the file [`CSUSHPINSA.csv`](https://github.com/AbirS2026/Double-LOcal-Unity/blob/main/CSUSHPINSA.csv) in the working directory.

Because Yahoo Finance and FRED data may be updated or revised over time, exact numerical outputs may differ slightly depending on the date on which the scripts are run.

---

## Related Applied Work

The historical bubble and date-stamping motivation is connected to the companion applied paper:

**Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance**  
Abir Sarkar and Martin T. Wells  
Paper: [arXiv:2604.12062](https://arxiv.org/pdf/2604.12062)

That paper focuses on robust date-stamping of bubble origination and collapse, while the present repository focuses on volatility-robust inference for autoregressive persistence and local explosiveness under nearly nonstationary stochastic volatility.

---

## Citation

If you use this repository, please cite the methodological paper:

```bibtex
@article{sarkar2025double,
  title={Double Local-to-Unity: Estimation under Nearly Nonstationary Volatility},
  author={Sarkar, Abir and Wells, Martin T},
  journal={arXiv preprint arXiv:2512.06823},
  year={2025}
}
```

Please also cite the companion applied paper when using the historical bubble and date-stamping motivation:

```bibtex
@article{sarkar2026there,
  title={Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance},
  author={Sarkar, Abir and Wells, Martin T},
  journal={arXiv preprint arXiv:2604.12062},
  year={2026}
}
```

The empirical data source is Yahoo Finance:

```bibtex
@online{yahooFinance,
  author  = {{Yahoo! Finance}},
  title   = {Yahoo historical data},
  year    = {2026},
  url     = {https://finance.yahoo.com}
}
```

---

## Notes

The code is organized so that each main figure and table can be reproduced from a corresponding script or notebook. File names are chosen to match the simulation and empirical components of the paper.

Yahoo Finance data are downloaded dynamically by the empirical scripts. Since adjusted prices may be revised over time, small differences from the paper's reported numbers may occur if the scripts are run at a later date.

For questions, please contact:

**Abir Sarkar**  
Cornell University  
`as4458@cornell.edu`
