- Andrea De Vita
- Enrico Lupi

-----------------------

# Bayesian SIR model with change points with application to the Omicron wave in Singapore

Final project for *Information Theory* course, Physics of Data degree at University of Padova.

## Abstract

In this work we reproduce the analysis conducted by J. Gu and G. Yin in their paper [Bayesian SIR model with change points with application to the Omicron wave in Singapore](https://www.nature.com/articles/s41598-022-25473-y). To analyze the impact of society or policy changes on the development of the Omicron wave, the stochastic susceptible-infected-removed (SIR) model with change points is proposed to accommodate the situations where the transmission rate and the removal rate may vary significantly at change points. Bayesian inference based on a Markov chain Monte Carlo algorithm is developed to estimate both the locations of change points as well as the transmission rate and removal rate within each stage, and experiments on simulated dataare carried out to test the effectiveness of the proposed method.

## Repository Structure

**Data** folder contains the simulated Markov chains for different scenarios and priors for &delta;.

**Presentation** folder contains the Jupyter notebook with example code and discussion of methodologies and results obtained.

**Scripts** folder contains the Python functions used in the project.
