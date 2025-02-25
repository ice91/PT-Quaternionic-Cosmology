# PT-Quaternionic-Cosmology-Colab

This repository provides code and data for testing the PT-symmetric quaternionic spacetime model in cosmology. It includes two main analyses:

1. **SN Ia Data Fitting:**  
   The Colab Notebook `SN_Ia_Fitting/pantheon_mcmc_colab.ipynb` demonstrates how to load the Pantheon SN Ia dataset, define the cosmological model (with an effective dark energy density parameter \(\Omega_{\mathrm{eff}}\)), and perform an MCMC analysis to obtain the distance-modulus vs. redshift plots and corner plots.

2. **GW170817 Analysis:**  
   The Colab Notebook `GW170817_Analysis/gw_speed_test_colab.ipynb` computes the relative deviation in gravitational-wave speed (\(\delta c_g/c\)) based on GW170817 and GRB170817A data.

## How to Use

You can run the notebooks directly in Google Colab:
1. Click on the respective notebook files.
2. Select **"Open in Colab"** (a button can be added via a link in the README).
3. Follow the instructions in the notebook to execute the code.

## Repository Structure

- **data/**: Contains the Pantheon SN Ia data file `lcparam_full_long_zhel.txt`.
- **SN_Ia_Fitting/**:
  - `pantheon_mcmc.py`: Original Python script for SN Ia data fitting (optional).
  - `pantheon_mcmc_colab.ipynb`: Colab Notebook for SN Ia data fitting.
- **GW170817_Analysis/**:
  - `gw_speed_test.py`: Original Python script for GW170817 analysis (optional).
  - `gw_speed_test_colab.ipynb`: Colab Notebook for GW170817 analysis.
- **docs/**: Contains additional documentation such as `Paper_Explanation.md`.

## Dependencies

The notebooks are set up to install any necessary packages automatically in the Colab environment. No additional setup is required.

## License

[ MIT License]
