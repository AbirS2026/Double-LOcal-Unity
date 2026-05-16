# Double Local-to-Unity: Estimation under Nearly Nonstationary Volatility

Replication code for the paper:

**Abir Sarkar and Martin T. Wells, “Double Local-to-Unity: Estimation under Nearly Nonstationary Volatility.”**  
Paper link: [https://arxiv.org/pdf/2512.06823](https://arxiv.org/pdf/2512.06823)

This repository contains the simulation and empirical code used to reproduce the main figures and tables in the paper. The paper develops volatility-robust moderate-deviation inference for autoregressive roots when both the mean dynamics and the stochastic-volatility process are highly persistent. The empirical applications study bubble diagnostics across historical exuberance episodes, global real-estate markets, commodity markets, and recent high-volatility comparison windows.

## Repository contents

### `Tables_1_2_KS_diagnostics.ipynb`

Generates the Monte Carlo Kolmogorov--Smirnov diagnostics reported in **Tables 1 and 2** of the paper.

The notebook simulates both:

- homoskedastic AR(1) innovations, and
- nearly nonstationary stochastic-volatility innovations.

It compares the simulated normalized OLS statistics with their theoretical limits:

- Gaussian limit in the nearly stationary regime,
- Cauchy limit in the mildly explosive regime.

Main outputs:

```text
table1_near_stationary_ks.csv
table2_mild_explosive_ks.csv
