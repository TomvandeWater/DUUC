""" Reference data from the fact sheet ATR72-600"""

""" Engines PW Canada PW127 M/N"""
P_TO = 2475  # Take off power [SHP]
P_TO_EO = 2750  # Take off power one engine [SHP]
P_max_con = 2500  # Power max continuous [SHP]
P_max_climb = 2192  # Power max climb [SHP]
P_max_cruise = 2132  # Power max cruise [SHP]

""" Blades Hamilton Standard 568F """
blade_diameter = 3.93  # blade diameter [m]
n_blades = 6  # number of blades [-]

""" Weights """
MTOW = 22800  # Maximum take-off weight (basic) [kg]
ZFW_max = 20800  # Max zero fuel weight (basic) [kg]
OEW = 13010  # Operational Empty weight (Tech spec.) [kg]
payload_max = 8550  # Max payload (at typical in-service OEW) [kg]
load_fuel = 5000  # Max fuel load [kg]

""" En-route performance"""
v_cruise = 510  # Max cruise speed (95% MTOW - ISA - Optimum FL) [km/h]
m_fuel = 762  # Fuel flow at cruise speed [kg/hr]
range = 758  # Range with max pax [NM]


