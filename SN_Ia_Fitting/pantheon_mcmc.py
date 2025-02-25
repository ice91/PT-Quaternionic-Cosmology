import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
import emcee
import corner

# Set a random seed for reproducibility
np.random.seed(42)

# -----------------------------------------------
# 1. Load the Pantheon SN Ia Data
# -----------------------------------------------
# Data file path (ensure the file "lcparam_full_long_zhel.txt" is placed in the "data" folder)
data_path = "../data/lcparam_full_long_zhel.txt"

# Read the file using whitespace as the delimiter, ignoring comment lines that start with '#'
sn_data = pd.read_csv(data_path, delim_whitespace=True, comment="#", header=None,
                      names=["name", "zcmb", "zhel", "dz", "mb", "dmb", "x1", "dx1",
                             "color", "dcolor", "3rdvar", "d3rdvar", "cov_m_s",
                             "cov_m_c", "cov_s_c", "set", "ra", "dec", "biascor"])
print("Successfully loaded SN Ia data. Data shape:", sn_data.shape)
print(sn_data.head())

# Use 'zhel' as redshift, 'mb' as observed magnitude, and 'dmb' as its error.
# Here, we define the observed distance modulus as μ_obs = mb - M, where M (absolute magnitude)
# is treated as a nuisance parameter.
z_data = sn_data["zhel"].values
mu_obs = sn_data["mb"].values  
mu_err = sn_data["dmb"].values

# Filter out data with z <= 0.01 to avoid large uncertainties at very low redshift.
mask = z_data > 0.01
z_data = z_data[mask]
mu_obs = mu_obs[mask]
mu_err = mu_err[mask]

# -----------------------------------------------
# 2. Define the Cosmological Model: Modified Friedmann Equation and Distance Modulus
# -----------------------------------------------
# We assume a flat universe (Omega_k = 0), with:
#   H(z) = H0 * sqrt(Omega_m * (1+z)^3 + Omega_eff)
# and distance modulus:
#   μ = 5 log10(d_L) + 25 + M_offset,
# where d_L is in Mpc, and M_offset accounts for the unknown absolute magnitude.
c = 299792.458  # Speed of light in km/s

def E_z(z, Omega_m, Omega_eff):
    return np.sqrt(Omega_m * (1 + z)**3 + Omega_eff)

def luminosity_distance(z, H0, Omega_m, Omega_eff):
    integral, _ = quad(lambda zp: 1.0 / E_z(zp, Omega_m, Omega_eff), 0, z)
    d_L = (1 + z) * (c / H0) * integral  # in Mpc
    return d_L

def mu_th(z, H0, Omega_m, Omega_eff, M_offset):
    # μ = 5 log10(d_L) + 25 + M_offset
    d_L = luminosity_distance(z, H0, Omega_m, Omega_eff)
    return 5 * np.log10(d_L) + 25 + M_offset

# -----------------------------------------------
# 3. Define the Log-Prior, Log-Likelihood, and Log-Probability Functions
# -----------------------------------------------
def log_prior(theta):
    # theta contains [H0, Omega_m, Omega_eff, M_offset]
    H0, Omega_m, Omega_eff, M_offset = theta
    if 50 < H0 < 90 and 0.1 < Omega_m < 0.5 and 0.5 < Omega_eff < 0.9 and -20 < M_offset < 0:
        return 0.0
    return -np.inf

def log_likelihood(theta, z, mu_obs, mu_err):
    H0, Omega_m, Omega_eff, M_offset = theta
    mu_model = np.array([mu_th(zi, H0, Omega_m, Omega_eff, M_offset) for zi in z])
    chi2 = np.sum(((mu_obs - mu_model) / mu_err)**2)
    return -0.5 * chi2

def log_probability(theta, z, mu_obs, mu_err):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, z, mu_obs, mu_err)

# -----------------------------------------------
# 4. Perform MCMC Parameter Fitting
# -----------------------------------------------
# Initial guess: H0 = 70 km/s/Mpc, Omega_m = 0.3, Omega_eff = 0.7, M_offset = -19.3
initial = np.array([70, 0.3, 0.7, -19.3])
ndim = len(initial)
nwalkers = 32
pos = initial + 1e-2 * np.random.randn(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(z_data, mu_obs, mu_err))

print("Starting MCMC fitting...")
nsteps = 5000  # Adjust the number of steps as needed
sampler.run_mcmc(pos, nsteps, progress=True)
print("MCMC fitting complete!")

samples = sampler.get_chain(discard=int(nsteps/2), flat=True)
print("Number of posterior samples:", samples.shape)

# -----------------------------------------------
# 5. Results Analysis and Visualization
# -----------------------------------------------
fig_corner = corner.corner(samples, labels=["$H_0$", "$\\Omega_m$", "$\\Omega_{eff}$", "$M_{offset}$"])
plt.show()

H0_mcmc, Omega_m_mcmc, Omega_eff_mcmc, M_offset_mcmc = np.median(samples, axis=0)
H0_err = np.percentile(samples[:,0], [16,84])
Omega_m_err = np.percentile(samples[:,1], [16,84])
Omega_eff_err = np.percentile(samples[:,2], [16,84])
M_offset_err = np.percentile(samples[:,3], [16,84])
print(f"H0 = {H0_mcmc:.2f} km/s/Mpc, 16/84% = {H0_err[0]:.2f}/{H0_err[1]:.2f}")
print(f"Omega_m = {Omega_m_mcmc:.3f}, 16/84% = {Omega_m_err[0]:.3f}/{Omega_m_err[1]:.3f}")
print(f"Omega_eff = {Omega_eff_mcmc:.3f}, 16/84% = {Omega_eff_err[0]:.3f}/{Omega_eff_err[1]:.3f}")
print(f"M_offset = {M_offset_mcmc:.2f}, 16/84% = {M_offset_err[0]:.2f}/{M_offset_err[1]:.2f}")

# -----------------------------------------------
# 6. Plot the Best-Fit Model vs. SN Ia Data
# -----------------------------------------------
def model_mu(z_array, H0, Omega_m, Omega_eff, M_offset):
    return np.array([mu_th(zi, H0, Omega_m, Omega_eff, M_offset) for zi in z_array])

z_plot = np.linspace(0.01, 1.8, 100)
mu_plot = model_mu(z_plot, H0_mcmc, Omega_m_mcmc, Omega_eff_mcmc, M_offset_mcmc)

plt.errorbar(z_data, mu_obs, yerr=mu_err, fmt='o', markersize=3, label='SN Ia Data', alpha=0.6)
plt.plot(z_plot, mu_plot, 'r-', lw=2, label='Best-fit Model')
plt.xlabel("Redshift $z$")
plt.ylabel("Distance Modulus $\\mu$")
plt.legend()
plt.title("SN Ia Distance Modulus vs. Redshift")
plt.show()

print("\nBest-fit results:")
print(f"H0 = {H0_mcmc:.2f} km/s/Mpc")
print(f"Omega_m = {Omega_m_mcmc:.3f}")
print(f"Omega_eff = {Omega_eff_mcmc:.3f}")
print(f"M_offset = {M_offset_mcmc:.2f}")
print("\nIf Omega_eff ≈ 0.7, then the dark energy density is consistent with our model's prediction (~10^-47 GeV^4).")
