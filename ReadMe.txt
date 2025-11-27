MEGB Simulations Project -------------------------------------------------------
This repository contains the code and scripts for model-based simulations using 
the MEGB (Mixed Effects Gradient Boosting) package.  
It provides implementations for multiple estimators, data generation scenarios, 
and evaluation of simulation results.  

The project is organized into two main model-based simulations as follows


Modelbased A ------------------------------------------------------------------
Contains scripts for running simulations under four different data generation 
scenarios using multiple estimators. Includes helper functions, estimation 
methods, and result aggregation scripts.

est.R 
– Central script to run all estimators

get_results.R 
– Aggregate results and calculate quality metrics

auxiliar
- contains helper functions and modules

estimators
- contents various estimation methods for simulations, e.g. BHF, EBP, GB, ...

data_generation_mod1.R – data_generation_mod4.R
- Generate synthetic datasets with different settings

Functions that implement the GBT, MERT and REEM
- gtb.R – Function for Gradient Tree Boosting 
- MEM_Boosting_VersionMERT.R
- MEM_Boosting_VersionREEM.R 


Modelbased B -------------------------------------------------------------------
Contains a simulation workflow with one data generation script, multiple 
estimator scripts, and result aggregation for evaluation.

run_simulations.R 
– Controls the full simulation workflow

results_sae.R 
– Aggregate results

data_generation.R 
– Generate synthetic dataset

QualityMeasure.R 
– Function to computes quality metrics

Functions that implement the MERT and REEM
- gtb.R – Function for Gradient Tree Boosting (Salditt et al.)
- MEM_Boosting_VersionMERT.R / MEM_Boosting_VersionREEM.R (Salditt et al.)
– Mixed-effects boosting


Project Structure --------------------------------------------------------------

├── MEGB_0.0.0.9000.tar.gz             R-package file (contains MEGB functions)
│
├── Modelbased A/                     Model-based simulations – Version A
│   ├── auxiliar/                      Helper functions and additional modules
│   │   ├── BHF_estimation.R           BHF estimator (Battese-Harter-Fuller)
│   │   ├── ebp.R                      Empirical Best Prediction (EBP)
│   │   ├── megb_package.R             Interface to the MEGB package
│   │   └── QualityMeasure.R           Calculation of quality measures
│   │
│   ├── estimators/                    Implementation of all estimators
│   │   ├── bhf.R                      BHF estimator (main implementation)
│   │   ├── BoostMERT.R                Boosting with MERT trees
│   │   ├── BoostREEM.R                Boosting with RE-EM trees
│   │   ├── ebp.R                      EBP estimator (estimator version)
│   │   ├── GB.R                       Basic Gradient Boosting
│   │   ├── MBOOST_T.R                 Model-based Boosting (tree variant)
│   │   ├── MBOOST.R                   General Model-based Boosting
│   │   ├── MEGB.R                     Mixed Effects Gradient Boosting
│   │   ├── MERF.R                     Mixed Effects Random Forest
│   │   ├── RF.R                       Random Forest (without random effects)
│   │   └── xgboost.R                  XGBoost implementation
│   │
│   ├── data_generation_mod1.R         Data generation – Scenario 1
│   ├── data_generation_mod2.R         Data generation – Scenario 2
│   ├── data_generation_mod3.R         Data generation – Scenario 3
│   ├── data_generation_mod4.R         Data generation – Scenario 4
│   ├── est.R                          Main estimation script (calls estimators)
│   ├── get_results.R                  Aggregation and calculation of results
│   ├── gtb.R                          Gradient Tree Boosting
│   ├── MEM_Boosting_VersionMERT.R     Mixed-effects boosting (MERT version)
│   └── MEM_Boosting_VersionREEM.R     Mixed-effects boosting (RE-EM version)
│
└── Modelbased B/                      Model-based simulations – Version B
    ├── data_generation.R              Data generation
    ├── gtb.R                          Gradient Tree Boosting (Salditt et al.)
    ├── MEM_Boosting_VersionMERT.R     Mixed-effects boosting (MERT version) (Salditt et al.)
    ├── MEM_Boosting_VersionREEM.R     Mixed-effects boosting (RE-EM version) (Salditt et al.)
    ├── QualityMeasure.R               Calculation of quality measures
    ├── results_sae.R                  Small Area Estimation results
    └── run_simulations.R              Controls the full simulation workflow