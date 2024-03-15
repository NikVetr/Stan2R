# Stan2R

Stan2R is a set of utilities I wrote to help convert [Stan](https://mc-stan.org/) models to R code. I did this for a few reasons -- sometimes I wanted to more easily perform prior or posterior predictive simulation or to estimate observation-wise posterior densities / masses, or I wanted to interactively step through the execution of a complicated model to figure out where something was going wrong. Occasionally I'd want to return to a fitted model and generate quantities for it, and did not want to use [generate_quantities](https://mc-stan.org/docs/cmdstan-guide/generate_quantities_config.html) in CmdStan (where one could much of the same functionality).

Also, collaborators are often unfamiliar with Stan, but were familiar with R. To help them understand a Stan model, I could trivially convert it to interactive R code for their ease of use.

Some of the functions have also come in handy elsewhere, eg in easily "flattening" multilevel models for fitting via methods besides MCMC (eg variational Bayes) as nodes in other approximation workflows (eg Empirical Bayes).

NOTE: do make sure to inspect any outputted R code for errors. I mostly wrote it with my own style preferences in mind, which may not match your own. I also did not exhaustively implement conversions for all Stan functions or distributions, or even for all features of Stan code (eg custom functions). Depending on the inputted model, you may have to modify converted scripts.

## Example

I wrote a short example (`\R\example.R`) using a fairly simple multilevel beta-binomial model.

In it, I perform:

- **Model Conversion**: Converts Stan model (`\models\beta_binomial.stan`) into executable R code.
- **Prior Predictive Simulation**: Generates samples from the model's prior predictive distribution, saving them to `\data\beta-binomial_simulated-data.json` and `\output\beta-binomial_simulated-params.json`.
- **Model Fitting**: Fits the model to a sample from this prior predictive distribution using CmdStanR.
- **Posterior Predictive Simulation**: Uses the fitted MCMC samples to simulate from the posterior predictive distribution.
- **Posterior Predictive Mass Estimation**: Uses the fitted MCMC samples to compute average fit mass for each observation.
- **Visual Assessment of Calibration**: Visually inspects model calibration by generating histograms of parameter quantiles in each of their marginal joint posterior distributions.
- **Visual Assessment of Model Fit**: Compares posterior predictive means against observed output values and unobserved variables (whose true values are known, as they were sampled from model prior). Visualizes these with scatterplots, alongside estimates of mean log posterior predictive masses and credible intervals.

This outputs a few quick figures, like

![Outcome Calibration](figures/outcome_calibration.png)

(showing quantiles for outcomes -- deviation from uniform(0,1) representing poor calibration)

![Parameter Fit](figures/parameter_goodness-of-fit.png)

(showing posterior means for different model parameters plotted against their true values)

## Getting Started

## Prerequisites

Before you begin, ensure you have R (version 4.0.0 or later) and CmdStan installed on your system. CmdStanR requires an installation of CmdStan to interface with Stan's C++ backend. Follow the instructions on the [CmdStanR GitHub page](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) for CmdStan installation guidelines.

## Installation

### 1. Install Required R Packages

Stan2R utilizes a set of R packages for data manipulation, Bayesian inference, and visualization. Install the required packages by running the following commands in your R environment:

```R
# Install CRAN packages
install.packages(c("data.table", "dplyr", "cmdstanr", "posterior", "extraDistr", "jsonlite", "parallel", "viridisLite"))

# Optionally, ensure cmdstanr is configured with CmdStan
cmdstanr::install_cmdstan()

## Installation

Clone the repository:

```bash
git clone https://github.com/NikVetr/Stan2R.git

## Usage

Please see `\R\example.R` for example usage -- you're probably best served by modifying that script for your own use. If things don't work, feel free to reach out and I'll do my best to help.

## Contributing

Contributions or extensions to this program are welcome! If you have suggestions for improvements or new features, please feel free to fork the repository and submit a pull request.