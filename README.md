# Job Scripts Repository

A collection of scripts for reproducing all displayed results.

## Adaptive_SPSA Scripts

### Purpose:

Reproduce the point estimates for simulated FHN data that gives the density plots (Fig. 4 and 7) .

### Scripts:

1.  **`bridge_4_005_LT.r`**
    -   What it does: For a given seed, simulates FHN data and maximizes the 4-step bridged Lie-Trotter partial likelihood using SPSA\
    -   Customise: Set seed number (Line 982)\
    -   Output: Numerical matrix of full convergence trajectory
2.  **`bridge_4_005_S.r`**
    -   What it does: For a given seed, simulates FHN data and maximizes the 4-step bridged Strang partial likelihood using SPSA\
    -   Customise: Set seed number (Line 1003)\
    -   Output: Numerical matrix of full convergence trajectory
3.  **`bridge_8_005_LT.r`**
    -   What it does: For a given seed, simulates FHN data and maximizes the 8-step bridged Lie-Trotter partial likelihood using SPSA\
    -   Customise: Set seed number (Line 982)\
    -   Output: Numerical matrix of full convergence trajectory
4.  **`bridge_8_005_S.r`**
    -   What it does: For a given seed, simulates FHN data and maximizes the 8-step bridged Strang partial likelihood using SPSA\
    -   Customise: Set seed number (Line 1004)\
    -   Output: Numerical matrix of full convergence trajectory

### Execution:

-   Each script is independent
-   Was run with seeds 1-100 for full results

## cSMC Package

### Purpose:

Contains the installable R package (`cSMC`) implementing a very basic cSMC; it is required for the MCMC and application script.

### Installation:

1.  Open `cSMC.Rproj` in RStudio
2.  Install the package in the `Build` interface

## MCMC Scripts

### Purpose:

Reproduce the posterior density plot (Fig. 3) for simulated FHN data.

### Dependencies

-   `mvnfast`, `MASS`, `mvtnorm`, `cSMC`

### Scripts:

1.  **`BPF_mcmc.R`**
    -   What it does: Produces result corresponds to the red line of Fig. 3
    -   Output: Numerical matrix containing the full BPF-PMMH trajectory
2.  **`cSMC_mcmc.R`**
    -   What it does: Produces result corresponds to the blue line of Fig. 3
    -   Output: Numerical matrix containing the full cSMC-PMMH trajectory

### Execution:

-   Both scripts are independent

## Cubic_SDE Script

### Pupose

Reproduces the convergence of likelihood result (Fig. 2) on the 1-d cubic SDE.

### Dependencies

-   `mvnfast`, `MASS`, `mvtnorm`, `doParallel`

## Application Script

### Pupose

Reproduces results related to the application on the rat (Fig. 8 and 9)

### Customise

-   Load the data that was originally in `.mat` files (Line 1-3)
### Dependencies

-   `mvnfast`, `MASS`, `mvtnorm`, `cSMC`
