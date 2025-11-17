# N2O_European_Coastal_Modelling
Spatial Analysis of in-situ N2O observations + RandomForest &amp; XGBoost to model coastal fluxes

# Project Overview

This repository contains all data, code, and model outputs associated with the study:

“Modeling Nitrous Oxide Fluxes in European Coastal Systems Using Machine Learning:
Analysing Regional Patterns, Drivers, and Climate Sensitivity.”

The project compiles ~8,500 in-situ N₂O flux measurements from European coastal waters (1993–2023), harmonizes them into a single dataset, and applies machine-learning models (Random Forest and XGBoost) to:
	•	Predict spatial distributions of coastal N₂O fluxes
	•	Identify dominant environmental drivers
	•	Benchmark against global ensemble estimates
	•	Conduct warming, oxygen decline, and nutrient perturbation scenarios
	•	Quantify prediction uncertainty
	•	Produce high-resolution (0.25°) gridded maps of coastal emissions

This repository ensures transparency and reproducibility for all analyses presented in the thesis and manuscript.


N2O_European_Coastal_Modelling/
│
├── data/
│   ├── Data_N2O_GuidedResearch.xlsx       # harmonised observational database
│   ├── *.nc                               # environmental predictor NetCDFs (WOA, Copernicus, GEBCO)
│   └── yang2020/                          # Yang et al. (2020) benchmark flux & dN2O files
│
├── scripts/
│   ├── 01_Overall_Map.R                   # Overall European Plots and data tidying and formatting
│   ├── 02_Zonal_Plot.R                    # Regional Summaries and Plots
│   ├── 03_Source_Sink.R                   # Identifying regions as sources or sinks of N2O
│   ├── 04_Modelling.R                     # Modelling script to predict N2O emissions using RF and XGBoost
│   ├── 05_Modelling_Updated.Rmd           # Updated and a cleaner version of the modelling script
|   ├── Bounding_Box.R                     # Create a visualisation of the bounding box for European coastal systems
| 
└── README.md


Environmental predictors (all long-term climatologies):
	•	World Ocean Atlas 2018/2023 – temperature, salinity, oxygen, nitrate, phosphate, density
	•	Copernicus Marine Service (CMEMS) – chlorophyll-a hindcast (1993–2025)
	•	GEBCO 2025 Grid – bathymetry
	•	Yang et al. (2020) global N₂O flux dataset – external benchmark

All .nc files are provided in the data/ directory unless restricted by license.

# Methods Overview

Machine Learning Models
	•	Random Forest (500 trees, tuned mtry)
	•	XGBoost (tuned with η=0.05, depth=6, subsample=0.8)
	•	Log-transformed flux target to reduce heteroscedasticity
	•	10-fold cross-validation for all models
	•	External benchmarking against Yang et al. (2020)
	•	SHAP values used for interpretability

Scenarios Modeled
	•	Warming: +2°C
	•	Deoxygenation: –20% O₂
	•	Nutrient enrichment: +20% NO₃⁻ & PO₄³⁻
	•	Combined multi-driver

Uncertainty Analysis
	•	100-member bootstrap Random Forest ensemble
	•	Standard deviation + relative uncertainty (%) maps


 # Key Findings (Summary)
	•	European coastal waters are consistent N₂O sources, especially the North Sea and Baltic Sea.
	•	Salinity, depth, chlorophyll, and temperature are dominant predictors of flux variability.
	•	Deoxygenation had the largest scenario impact (+10%).
	•	Combined climate and nutrient forcing increases emissions by ~13%.
	•	Prediction uncertainty is lowest in well-sampled basins and highest in data-poor regions such as the Eastern Mediterranean and Black Sea.