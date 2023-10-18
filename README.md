# combined-quantile-forecast

The source code for simulations and real data analysis of $\text{PM}_{2.5}$ data in South Korea, included in the paper 'Combined quantile forecasting for high-dimensional non-Gaussian data'.

## Overview

- Codes
  - `functions.R` is a code for functions used in analysis, including six forecasting methods. 
  - `Simulation.R` is a code that implements a simulation study and derives the results in Table 1 and 2.
  - `Simulation_extend.R` is a code of simulation study for our extended method that derives the results in Table 3.
  - `PM_proposed.R` contains the process of real data anaylsis of PM2.5 data using the proposed method. Results are included in Table 5 and 6.
  - `PM_compete.R` contains the process of real data anaylsis of PM2.5 data using competing methods. Results are included in Table 5 and 6.
  
- Data
  - `PM.csv` contains $\text{PM}_{2.5}$ concentration observed for 20784 hours in 308 stations in South Korea. The missing values are already imputed by the EM algorithm mentioned in the paper.
  - `station_loc.csv` contains information about the location of 308 stations. 
  - `r.vec_list.rds` contains a single list of the estimated number of quantile factors for each period. It can be used to skip the estimation process which is included in the `PM_proposed.R` code. 
 
