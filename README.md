# Stochastic-Modeling-NeutronFlux
# Stochastic Time Series Analysis of Galactic Cosmic Ray Intensity
### Radiation Hazard Forecasting for the Alcomsat-1 Satellite Mission

---

## Table of Contents
1. [Project Overview](#project-overview)
2. [Scientific Background](#scientific-background)
3. [Repository Structure](#repository-structure)
4. [Datasets](#datasets)
5. [Methodology](#methodology)
6. [Models](#models)
7. [Key Results](#key-results)
8. [SEU Hazard Index](#seu-hazard-index)
9. [Requirements](#requirements)
10. [How to Run](#how-to-run)
11. [Authors & Context](#authors--context)

---

## Project Overview

This project applies stochastic time series analysis to **Galactic Cosmic Ray (GCR) intensity**, measured as monthly neutron flux deviation at the Oulu Neutron Monitor Station (Finland), covering the period **March 1964 to August 2025** (737 months).

The primary objective is to forecast the space radiation environment relevant to the **Alcomsat-1** Algerian geostationary communications satellite. High-energy GCR particles pose a direct operational risk to satellite electronics through **Single Event Upsets (SEUs)** — commonly known as *bit flips* — where a particle strike deposits enough charge to flip a binary state in a microprocessor or memory chip. Such events can corrupt telemetry data, cause software crashes, or trigger unintended commands affecting satellite orientation.

The project delivers:
- A systematic comparison of **6 time series and statistical models**
- A diagnosis of why standard ARIMA differencing fails for physically-driven cyclical data
- An **operational SEU Hazard Index** $H(t)$ translating flux forecasts into actionable risk probabilities for mission controllers

---

## Scientific Background

### Galactic Cosmic Rays (GCRs)
GCRs are high-energy particles — primarily protons — originating from outside the solar system. Their intensity at Earth is not constant: it is modulated by the Sun's magnetic field over the 11-year solar cycle (Schwabe cycle).

### The Inverse Relationship (Solar Modulation)
The fundamental physical mechanism governing this dataset:

$$\text{Solar Activity (SSN)} \uparrow \;\Rightarrow\; \text{Heliospheric Magnetic Field} \uparrow \;\Rightarrow\; \text{GCR Deflection} \uparrow \;\Rightarrow\; \text{Neutron Flux} \downarrow$$

- **Solar Maximum**: Strong, turbulent magnetic field acts as a shield. Neutron flux is at minimum.
- **Solar Minimum**: Weakened magnetic shield allows more GCRs to penetrate. Neutron flux peaks.

This produces the characteristic slow, smooth, ~132-month oscillation visible in the data, with a measured Pearson correlation of **r = −0.81** between Sunspot Number (SSN) and neutron flux.

### Why Neutron Flux?
GCRs are difficult to measure directly from the ground. When a GCR strikes Earth's atmosphere it produces a cascade of secondary particles including neutrons. Ground-based neutron monitors provide a highly accurate and continuous proxy for the space GCR environment — far longer records than direct space measurements allow.

### Single Event Upsets (SEUs)
When a GCR or secondary particle strikes a microelectronic component, it deposits a localised charge. If sufficient, this flips a binary state (0 → 1 or 1 → 0). For Alcomsat-1:
- **High neutron flux period** (solar minimum) → elevated SEU rate → higher risk of bit flips in onboard computers
- **Low neutron flux period** (solar maximum) → suppressed SEU rate → safer window for critical operations

---

## Repository Structure

```
neutron-flux-timeseries/
│
├── Neutron_Flux_TimeSeries.ipynb   # Main analysis notebook (R kernel)
│
├── data/
│   ├── oulu_raw2.txt               # Oulu neutron monitor raw data
│   └── SN_m_tot_V2.0.csv          # SIDC international sunspot number
│
└── README.md                       # This file
```

---

## Datasets

### 1. Oulu Neutron Monitor — Neutron Flux
| Field | Value |
|-------|-------|
| **Source** | University of Oulu, Finland |
| **URL** | https://cosmicrays.oulu.fi |
| **Variable** | Neutron flux as % deviation from long-term mean |
| **Period** | 1964 – 2025 |
| **Frequency** | Monthly averages |
| **Format** | Semicolon-delimited text; columns: DateTime, NeutronFlux |

The Oulu station is one of the longest-running and most-cited neutron monitor stations in the world, providing a continuous, calibrated record spanning more than 6 complete solar cycles.

### 2. SIDC International Sunspot Number v2.0 — SSN
| Field | Value |
|-------|-------|
| **Source** | Solar Influences Data Analysis Center (SIDC), Royal Observatory of Belgium |
| **URL** | https://www.sidc.be/SILSO/datafiles |
| **File** | SN_m_tot_V2.0.csv |
| **Variable** | Monthly international sunspot number |
| **Period** | 1749 – present |
| **Frequency** | Monthly |
| **Note** | Version 2.0 applies a full historical recalibration relative to v1.0 |

Both datasets are merged on Year–Month. The analysis uses **March 1964 – December 2015** for training and **January 2016 – August 2025** (116 months) as a withheld test set.

---

## Methodology

### Train / Test Split
| Split | Period | Months |
|-------|--------|--------|
| Training | March 1964 – December 2015 | 621 |
| Test (withheld) | January 2016 – August 2025 | 116 |

The test set was withheld from all model training and used exclusively for out-of-sample evaluation.

### Framework I — Box-Jenkins ARIMA

The standard three-phase Box-Jenkins methodology:

**Phase I — Identification**
- Augmented Dickey-Fuller (ADF) test for stationarity
- Raw series: p = 0.73 (non-stationary)
- First-differenced series: p < 0.01 (stationary)
- ACF/PACF analysis to identify p, q orders
- Periodogram analysis confirming dominant period of 132 months (11-year solar cycle)

**Phase II — Estimation & Diagnostics**
- Maximum likelihood estimation of ARIMA parameters
- Ljung-Box test on residuals (H₀: residuals are white noise)
- QQ plots and ACF/PACF of residuals

**Phase III — Forecasting**
- Multi-step out-of-sample forecasting (116 months)
- Rolling 1-step-ahead walk-forward validation (24 months)
- Accuracy metrics: MAE, RMSE, Pearson correlation

### Framework II — Generalized Additive Model (GAM)

Introduced to address a fundamental limitation of the ARIMA family: the inability to preserve the **nonlinear SSN–flux relationship in levels** when differencing is applied.

**The differencing problem — quantified:**
Because flux depends on SSN, applying d=1 automatically converts the model from estimating the level relationship (r = −0.81) to the change relationship (r = −0.19):

$$\Delta\text{flux}(t) = \beta \cdot \Delta\text{SSN}(t) + \Delta\varepsilon(t)$$

This destroys the physically meaningful signal. The GAM avoids differencing entirely by fitting the relationship in levels with a nonlinear smooth function.

**Stationarity diagnosis:**
OLS regression of flux on SSN followed by ADF test on residuals gives p = 0.023 — confirming that **SSN is the source of non-stationarity**, not a stochastic unit root. This validates the d=0 approach.

---

## Models

| Model | Framework | Key Parameters | Notes |
|-------|-----------|---------------|-------|
| ARIMA(0,1,2) | Box-Jenkins | d=1, q=2 | Baseline; flat long-horizon forecast |
| SARIMA(0,1,2)(1,0,1)[132] | Box-Jenkins | Seasonal period=132 | Encodes 11-year cycle |
| ARIMAX d=1 | Box-Jenkins | d=1, xreg=SSN | Differencing destroys SSN signal |
| ARIMAX d=0 | Box-Jenkins | d=0, xreg=SSN | AR(1)=0.977 absorbs SSN signal |
| Rolling 1-Step ARIMA | Box-Jenkins | Walk-forward | Best short-term operational performance |
| GAM v2 | GAM | s(SSN, k=10) | Nonlinear SSN spline; no time trend |
| **GAM v3** | **GAM** | **s(SSN) + sin/cos cycle phase** | **Best model overall** |

### GAM v3 Specification

$$\text{flux}(t) = \alpha + f(\text{SSN}(t)) + \gamma_1 \sin\!\left(\frac{2\pi t}{132}\right) + \gamma_2 \cos\!\left(\frac{2\pi t}{132}\right) + \varepsilon(t)$$

where $f(\cdot)$ is a cubic regression spline (edf = 4.1, confirming nonlinearity). The harmonic terms encode solar cycle phase in a form that extrapolates periodically and correctly outside the training range — unlike polynomial time splines, which diverge.

---

## Key Results

### Out-of-Sample Performance (Jan 2016 – Aug 2025)

| Model | MAE | RMSE | Corr |
|-------|-----|------|------|
| ARIMA(0,1,2) | 6.844 | 7.362 | 0.010 |
| SARIMA(0,1,2)(1,0,1)[132] | 6.676 | 7.156 | 0.684 |
| ARIMAX d=1 (differenced) | 6.797 | 7.291 | −0.495 |
| ARIMAX d=0 (levels) | 6.685 | 7.117 | −0.495 |
| GAM v2 (SSN only) | 3.391 | 4.174 | 0.748 |
| **GAM v3 (SSN + cycle phase)** | **2.837** | **3.490** | **0.867** |

**GAM v3 reduces MAE by 58% relative to the best ARIMA-family model** (SARIMA) and achieves a directional correlation of 0.867, correctly tracking the full rise and fall of Solar Cycle 25 across the 116-month test period.

### Key Scientific Finding
The ARIMA models' failure is not due to model misspecification in the classical sense — all pass the Ljung-Box white noise test. The failure is **physical**: differencing converts a level-based relationship (r = −0.81) into a change-based relationship (r = −0.19), destroying the solar modulation signal that makes SSN useful as a predictor. The GAM framework avoids this by working entirely in levels.

---

## SEU Hazard Index

The GAM v3 forecast is converted into an operational **Single Event Upset Hazard Index**:

$$H(t) = \Phi\!\left(\frac{\hat{Y}(t) - \mu_Y}{\sigma_Y}\right)$$

where $\hat{Y}(t)$ is the GAM v3 flux forecast, $\mu_Y$ and $\sigma_Y$ are the historical mean and standard deviation, and $\Phi(\cdot)$ is the standard normal CDF. This maps flux forecasts to a probability in $[0, 1]$.

### Operational Traffic Light

| $H(t)$ | Risk Level | Recommended Action |
|--------|------------|-------------------|
| ≥ 0.75 | 🔴 HIGH | Defer software updates; avoid critical memory writes |
| 0.40 – 0.75 | 🟡 MODERATE | Increase error-correction polling frequency |
| < 0.40 | 🟢 LOW | Normal operations; schedule critical maneuvers |

---

## Requirements

The notebook runs on an **R kernel** (Jupyter or RStudio). Required packages:

```r
install.packages(c(
  "tseries",      # ADF test, time series utilities
  "forecast",     # auto.arima, Arima, forecast objects
  "mgcv",         # GAM fitting (gam, s())
  "lubridate"     # date handling
))
```

Tested on R 4.3+. No Python dependencies.

---

## How to Run

1. Clone or download this repository
2. Ensure R 4.3+ is installed with the packages listed above
3. Open `Neutron_Flux_TimeSeries.ipynb` in Jupyter (with R kernel) or convert to an RMarkdown file for RStudio
4. Run cells sequentially from top to bottom — each cell depends on objects created by previous cells
5. The data is loaded automatically from the GitHub raw URLs at the top of the notebook; no local data files are required to reproduce results
6. Identical data files are archived in the submitted Google Drive folder for offline verification

> **Note:** The SARIMAX cell with `stepwise=FALSE, approximation=FALSE` can take 10–30 minutes due to the large seasonal period (132). Set `stepwise=TRUE, approximation=TRUE` for a fast equivalent result.

---

## Authors & Context

Developed as a course project in **Stochastic Processes / Time Series Analysis**.

**Satellite context:** Alcomsat-1 is Algeria's first geostationary communications satellite, launched in December 2017 and operated by Algérie Télécom Satellite (ATS). It operates in geostationary orbit at 24.8°E, where GCR-induced SEUs represent a non-trivial operational risk over the satellite's designed 15-year lifespan.

**Data acknowledgements:**
- Oulu Neutron Monitor data: University of Oulu Cosmic Ray Station, Finland
- Sunspot Number: SILSO World Data Center, Royal Observatory of Belgium, Brussels
