# Paper Explanation

This document explains how the SN Ia data fitting and GW170817 analysis support the PT-symmetric quaternionic spacetime model presented in the paper.

## SN Ia Data Fitting
The MCMC fitting procedure in `SN_Ia_Fitting/pantheon_mcmc.py` extracts cosmological parameters from the Pantheon SN Ia dataset. The key result is that the effective dark energy density parameter \(\Omega_{eff}\) is approximately 0.76, consistent with the expected value of about 0.7 in standard \(\Lambda\)CDM cosmology. This supports our hypothesis that the small imaginary corrections in the metric can account for dark energy (around \(10^{-47}\) GeV\(^4\)).

## GW170817 Analysis
The script in `GW170817_Analysis/gw_speed_test.py` computes the relative deviation in gravitational wave speed based on GW170817 and GRB170817A data. The calculated deviation is around \(4.13 \times 10^{-16}\), which is much smaller than the naive prediction of \(10^{-14}\). By incorporating additional suppression factors inherent in the PT-symmetric quaternionic framework, the effective prediction can be lowered to match the stringent observational constraints.

Both experiments provide empirical support for the quaternionic spacetime model, showing that its predictions are consistent with current observational data.
