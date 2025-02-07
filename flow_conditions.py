from analysis_modules.aerodynamic import speed_of_sound
from analysis_modules.ISA import air_density_isa
import numpy as np

""" Flow parameters """
rho = 0.125  # air density at sea level [kg/m^3]
airspeed = 148  # [m/s]
u_inf = 180  # free stream flow velocity [m/s]
altitude = 7000  # flying altitude in [m]

""" Calculated flow parameters """
Mach = u_inf / (speed_of_sound(altitude))  # free stream mach number [-]

""" Flight conditions"""
""" dive speed is now a set number based on Stavreva, could potentially be 
a variable where VD => 1.25 V"""
V_d = 193.2  # dive speed [knots]

""" alpha must be between -5 and 15 otherwise there will be issues with pylon 
forces """
alpha = 5  # angle of attack [deg]
delta_e = 0  # elevator deflection angle [deg]
delta_r = 0  # rudder deflection angle [deg]


""" print section"""
print("\n----- Flow parameters -----")
print("Airspeed =", u_inf, "[m/s]")
print("Altitude =", altitude, "[m]")
print("Mach =", np.round(Mach), '[-]')
print("\n ----- Angles -----")
print("Angle of Attack:", alpha, "[deg]")
print("Elevator deflection", delta_e, "[deg]")
print("Rudder deflection", delta_r, "[deg]")
print("\n----- ISA calculator -----")
print(f'Air density =', np.round(air_density_isa(altitude)[0], 3), "[kg/m3]")
print(f'Air temperature =', np.round(air_density_isa(altitude)[1], 1) , "[K]")