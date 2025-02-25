import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. Define Observational Data for GW170817 and GRB170817A
# ------------------------------------------------------------------------------
# GW170817 (binary neutron star merger) was detected by LIGO/Virgo.
# GRB 170817A (associated gamma-ray burst) was observed approximately 1.7 seconds later.
# The estimated distance to the event is about 40 Mpc.
t_GW = 0.0         # Reference time for GW170817 (in seconds)
t_GRB = 1.7        # Arrival time of GRB 170817A (in seconds)
delta_t = t_GRB - t_GW  # Time delay between GRB and GW (in seconds)

distance_Mpc = 40.0           # Distance in Megaparsecs (Mpc)
Mpc_to_km = 3.086e19          # Conversion: 1 Mpc = 3.086e19 km
distance_km = distance_Mpc * Mpc_to_km

# Speed of light in km/s
c = 3.0e5  # km/s

# ------------------------------------------------------------------------------
# 2. Compute the Relative Deviation in Gravitational Wave Speed
# ------------------------------------------------------------------------------
# The light travel time is:
light_travel_time = distance_km / c  # seconds
# The relative deviation is given by:
relative_deviation = delta_t / light_travel_time

print("Calculated relative deviation in gravitational wave speed (δc_g / c):")
print(f"{relative_deviation:.2e}")

# ------------------------------------------------------------------------------
# 3. Visualize the Result
# ------------------------------------------------------------------------------
plt.figure(figsize=(6,4))
plt.bar(["δc_g / c"], [relative_deviation], color="blue")
# Draw a horizontal dashed line at 1e-14 as the naive theoretical prediction
plt.axhline(y=1e-14, color='red', linestyle='--', label="Predicted ~1e-14")
plt.yscale("log")
plt.ylabel("Relative Deviation (δc_g / c)")
plt.title("GW170817: Gravitational Wave Speed Deviation Test")
plt.legend()
plt.show()
