# Library imports
import numpy as np
import matplotlib.pyplot as plt

#####################################################################
# Initialization
#####################################################################
# Advance ratio
mu_arr = [0.0, 0.1, 0.25, 0.5]

# Induced inflow ratio
lambda_i = 0.09781  # averaged out over the radial direction

#####################################################################
#  Thrust ratio computations
#####################################################################


def thrust_ratio_calc(elevation, mu, lambda_i):
    correction_factor = 1 + (mu / lambda_i) ** 2
    thrust_ratio = (
        1
        + 0.0199 / (elevation ** 2 * correction_factor)
        + 0.0944 / (0.6572 + (elevation ** 2 * correction_factor))
    )
    return thrust_ratio


#####################################################################
# Calculations for discretized variables
#####################################################################
dis_num = 100  # Discretization

elevation_arr = np.linspace(0.2, 3, dis_num)
thrust_ratio_for_mu = []

for mu in mu_arr:
    thrust_ratio_for_mu.append(
        [
            thrust_ratio_calc(elevation, mu, lambda_i)
            for elevation in elevation_arr
        ]
    )

#####################################################################
# Post-processing
#####################################################################
# fig, [ax1, ax2] = plt.subplots(1, 2)
# fig.suptitle("Thrust vs Elevation")

# ax1.set_title("Without modifications")
# ax1.plot(elevation_arr, thrust_ratio_arr)
# ax1.set_xlabel("elevation (z/R)")
# ax1.set_ylabel("Thrust (IGE/OGE)")

# ax2.set_title("With modifications")
# ax2.plot(elevation_arr, thrust_ratio_arr_with_mod)
# ax2.set_xlabel("elevation (z/R)")
# ax2.set_ylabel("Thrust (IGE/OGE)")

# fig.tight_layout()
# plt.show()

plt.figure()
plt.title("Variation in thrust for forward flight", fontweight="bold")

for i in range(len(mu_arr)):
    plt.plot(
        elevation_arr, thrust_ratio_for_mu[i], label=f"$ \mu$ = {mu_arr[i]}"
    )

plt.legend()
plt.grid()
plt.xlabel("elevation (z/R)")
plt.ylabel("Thrust (IGE/OGE)")
plt.xlim((0, elevation_arr[-1]))
plt.ylim(bottom=0.9)
plt.show()

