# Historic vs Implied Volatility

**Report (PDF):** [[Zenodo DOI](https://doi.org/10.5281/zenodo.19674612)]  
**Code Repository:** [[GitHub link](https://github.com/jacekschmidt/Historial_vs_Implied_Volatility/edit/main/README.md)]

## Introduction

In every quantitative model of an asset price, two components determine its dynamics: a deterministic drift term and a stochastic volatility term. The drift, or trend, represents the expected average direction of the asset; for example, this may correspond to a long-term upward trajectory or to cyclical behaviour, depending on the data and the chosen specification. Because the drift is highly sensitive to the sample period and modelling choices, its estimation is often subjective. The stochastic component, volatility, measures the random deviations from this trend. It governs the magnitude of fluctuations in asset prices and is directly linked to the risk borne by investors. In share and derivatives markets, volatility plays a central role: higher volatility implies a wider distribution of possible outcomes for an asset, which translates into greater uncertainty in the valuation of options. As such, volatility estimation is the cornerstone of both risk management and derivative pricing. Beyond derivative markets, volatility is also critical in applications such as portfolio allocation, Value-at-Risk (VaR) calculations, and stress testing.

To study volatility, a model of asset price dynamics must first be specified. The choice of model is crucial because the decomposition into drift and volatility is not unique: a more refined drift specification can capture more systematic movements of the asset, leaving a smaller residual volatility term to represent only the truly random fluctuations. Conversely, if the drift is oversimplified, volatility will appear larger because it must compensate for unmodelled systematic behaviour. In this project, we adopt the geometric Brownian motion (GBM) model. GBM is widely used because it captures several stylised facts of financial markets: asset prices remain strictly positive, volatility scales with both time and price, and log returns are approximately normally distributed over short horizons. At the same time, GBM has important shortcomings. It does not capture volatility clustering nor sudden jumps, but it remains a natural starting point due to its analytical tractability and its role in option pricing.

Two main approaches to estimating volatility are considered: historical volatility and implied volatility. Historical volatility is backward-looking, relying on past price observations under the real-world probability measure. Implied volatility is forward-looking, inferred from option prices under the risk-neutral measure, and incorporates both market expectations and a volatility risk premium. The aim of this project is to examine how these two methods are formulated, evaluate their numerical implementation, and assess their reliability in practical applications.

---

## Project Structure

This project implements a numerical framework in C++ for estimating and comparing historical and implied volatility.

The codebase includes:

* Numerical integration routines (Simpson’s rule)
* Option pricing via transformed integrals
* Historical volatility estimation from time series data
* Root-finding algorithms for implied volatility
* Supporting mathematical utilities via custom vector classes

Header files used in this project:

mvector.h
cmath
fstream
vector
algorithm

---

## Mathematical Formulation

The asset price is modelled using geometric Brownian motion.

Historical volatility is computed from log returns:

σ ≈ sqrt(252 × variance of log returns)

Implied volatility is defined implicitly through the equation:

Option Price (σ) = Market Price

This requires solving a nonlinear equation, typically via iterative numerical methods.

---

## Algorithms Implemented

The project implements several numerical methods.

Simpson’s Rule is used for numerical integration of option pricing integrals. This method provides fourth-order accuracy and is efficient for smooth integrands.

The Newton-Secant method is used to compute implied volatility. This iterative method avoids explicit derivative computation while retaining fast convergence.

A variance estimator computes historical volatility using log returns derived from asset price data.

---

## Code Overview

The implementation is organised into several key components.

Payoff functions define the structure of European call and put options. These determine the terminal value of the option based on the transformed integration variable.

An abstract Integral class provides a framework for defining integrands and residual functions. Derived classes such as European_Call and European_Put implement the pricing integrals.

The integrate function applies Simpson’s rule to approximate definite integrals. This is used to compute option prices numerically.

The Variance_Estimator function calculates annualised volatility from historical price data by computing log returns and their variance.

The NewtonSecant function solves nonlinear equations to determine implied volatility. It iteratively updates the volatility estimate until the pricing error falls below a specified tolerance.

---

## Running the Code

To compile the program:

g++ -std=c++11 main.cpp -o volatility

To run the program:

./volatility

---

## Numerical Experiments

The project includes several computational experiments.

Historical volatility estimation is performed using real price data read from a file. The resulting volatility reflects past market behaviour.

Option pricing is evaluated using numerical integration, allowing comparison between computed values and known analytical results.

Implied volatility is computed by matching model prices to observed market prices using the Newton-Secant method.

Error analysis investigates how numerical integration accuracy depends on step size and domain truncation.

---

## Results and Observations

Historical volatility is simple to compute and requires only past data, but it is inherently backward-looking and may not capture current market expectations.

Implied volatility reflects forward-looking information embedded in option prices. It is more relevant for pricing and risk management but requires numerical computation and depends on model assumptions.

The Newton-Secant method demonstrates fast convergence when initial guesses are reasonable, though it may struggle if the function is poorly conditioned.

Simpson’s rule provides accurate integration for smooth functions, but computational cost increases with the number of intervals.

The geometric Brownian motion model simplifies analysis but introduces modelling error in markets exhibiting volatility clustering or jumps.

---

## Conclusion

This project demonstrates how historical and implied volatility can be computed within a unified numerical framework. Historical volatility offers a straightforward estimate based on past data, while implied volatility incorporates market expectations and is more relevant for derivative pricing. The implementation highlights the role of numerical integration and root-finding in financial modelling, as well as the trade-offs between simplicity, accuracy, and computational cost.

---

## Future Work

Potential extensions include implementing stochastic volatility models, improving root-finding robustness, and applying the framework to real-world datasets with multiple maturities and strike prices. Incorporating more advanced numerical methods and alternative pricing models would further enhance accuracy and applicability.

---

## Author

Jacek Schmidt
jacek14schmidt@gmail.com
