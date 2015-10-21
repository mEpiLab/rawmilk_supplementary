# Supplementary material for Rawmilk paper

This repository contains the necessary material to reproduce the paper
*Estimating bacterial pathogen concentrations in New Zealand bulk tank milk*.

## Supplementary text

The supplementary_material folder contains a PDF file of supplementary details which includes the
full formulation of the MCMC scheme for fitting total bacterial count data.

## Data

The data necessary for the total bacterial count model and pathogen survey simulation
are located in the data folder.

## Total Bacterial Counts versus other covariates

The generalised linear mixed model for fitting TBC to other covariates, and computation
of $R^2$ values is available in the tbc_vs_other_variables.R script in the tbc_vs_other_variables
folder.

## Total Bacterial Counts

The model for total bacterial counts are fitted to the data using the MCMC algorithm
described in tbc_model.R in the tbc_model folder.

Diagnostics may be produced using tbc_diagnostics.R script in the same folder.

## Pathogen survey simulation

The survey simulation, including predicted contentrations per litre is found in the
survey.R script in the survey_model folder.

## Figures

The figures in the paper may be reproduced, once the above models are run and output produced,
using figures.R located in the figures sub-folder.
