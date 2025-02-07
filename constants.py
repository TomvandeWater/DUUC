"""Constants"""

g = 9.81  # gravitational constant [m/s^2]
P0 = 101325  # Sea level standard atmospheric pressure (Pa)
T0 = 288.15  # Sea level standard temperature (K)
L = 0.0065  # Temperature lapse rate (K/m)
R = 287.05  # Specific gas constant for dry air (J/(kg·K))
mu_air = 1.7894e-5  # dynamic viscosity at sea level [m^2/s]
gamma_air = 1.4  # specific heat ratio of air


""" Constants in design"""
""" Pylon factor is based on estimated and calculated pylon weight 115/90 """
K_pylon = 1.30  # factor to compensate for underestimating pylon weight

""" value to determine weight of elevators empirically """
K_weight_h_elevator = 5  # in [kg/m2]
K_weight_v_elevator = 5  # in [kg/m2]
