""" Flow parameters """
rho = 0.125  # air density at sea level [kg/m^3]
airspeed = 148  # [m/s]
a = 1.2  # speed of sound at sea level [m/s]

""" Flight conditions"""
""" dive speed is now a set number based on Stavreva, could potentially be 
a variable where VD => 1.25 V"""
V_d = 193.2  # dive speed [knots]

""" alpha must be between -5 and 15 otherwise there will be issues with pylon 
forces """
alpha = 5  # angle of attack [deg]
Mach = 0.4  # free stream mach number [-]
Re = 100000  # Reynolds number [-]
u_inf = 180  # free stream flow velocity [m/s]
