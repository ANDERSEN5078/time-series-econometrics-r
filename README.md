# 📊 Time Series Econometrics in R

## 🎯 Project Overview

This project provides a complete econometric analysis of time series data using R.

It covers both theoretical concepts and practical implementation, including:

* Stationarity analysis
* Simulation of stochastic processes
* Spurious regression
* Cointegration
* Error Correction Models (ECM)

---

## 📦 Dataset

The data used in this project comes from the **PoEdata** package.

---

## ⚙️ Technologies & Libraries

* R
* tseries
* dynlm
* forecast
* lmtest
* sandwich
* broom
* knitr
* car

---

## 🔬 Methodology

### 1. Stationarity Analysis

* Visualization of time series
* Autocorrelation function (ACF)
* Augmented Dickey-Fuller (ADF) test

### 2. Simulation of Processes

* AR(1) processes
* Random Walk
* Random Walk with drift and trend

### 3. Spurious Regression

* Illustration of misleading regression results between non-stationary variables

### 4. Cointegration Analysis

* Long-run relationship estimation
* Residual stationarity testing
* Phillips-Ouliaris test

### 5. Error Correction Model (ECM)

* Short-run dynamics modeling
* Speed of adjustment estimation
* Nonlinear estimation using NLS

---

## 🚀 How to Run the Project

```r
install.packages("remotes")
remotes::install_github("ccolonescu/PoEdata")

source("analysis.R")
```

---

## 📈 Key Insights

* Non-stationary time series can lead to spurious regression results
* Differencing helps achieve stationarity
* Cointegration reveals long-term equilibrium relationships
* ECM models capture both short-term dynamics and long-term adjustments

---

## 📁 Project Structure

```
.
├── scripts/
│   └── analysis.R
├── README.md
```

---

## 👨‍💻 Author

**Lysias Andersen**

---

## 🌍 Perspective

This project is part of a broader journey toward mastering:

* Econometrics
* Data Science
* Artificial Intelligence

---

## ⭐ If you like this project

Feel free to star the repository and share it!
