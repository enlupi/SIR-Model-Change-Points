- Andrea De Vita
- Enrico Lupi

-----------------------

# Bayesian SIR model with change points with application to the Omicron wave in Singapore

Final project for *Information Theory* course, Physics of Data degree at University of Padova.

## Abstract

In this work we reproduce the analysis conducted by J. Gu and G. Yin in their paper [Bayesian SIR model with change points with application to the Omicron wave in Singapore](https://www.nature.com/articles/s41598-022-25473-y#:~:text=To%20analyze%20the%20impact%20of,vary%20significantly%20at%20change%20points.). 

## Repository Structure

The aforementioned files are available in the **References** folder, while the analyzed INFN dataset is contained in **datasets** as Excel files.

As the name suggests, the **SingleRun** folder contains the results of the analysis conducted on a single run of data taking, i.e. with only one setting for the cavity frequency. Inside, **prepData.ipynb** explains the procedure to take and prepare the data from the raw Excel files for the following analysis,  **fitFunc.ipynb** illustrates the functions for the background and signal fit while **Statistics.ipynb** describes the actual statistical procedure; lastly, **SingleRunAnalysis.ipynb** shows the obtained results.

The **MultipleRunAnalysis** folder cobtains the same elements and analysis as **SingleRun**, but adapts them to work simultaneously on the whole dataset with different cavity frequency settings.
