""" Propulsive empennage set up"""
duct_diameter = 3.6  # duct diameter in [m]
AR_duct = 2  # aspect ratio of the duct [-]
duct_chord = duct_diameter / AR_duct  # duct chord in [m]
ai_duct = 0  # installation angle duct [deg]

pylon_length = 1  # pylon length in [m]
pylon_chord = 1  # pylon chord length [m]
cant_angle = 30  # cant angle of the pylon wrt to the fuselage [deg]

nacelle_length = 1.3  # nacelle length in [m]
nacelle_diameter = 1.0  # nacelle diameter in [m]

support_length = duct_diameter  # support length in [m]
support_chord = 0.5  # support chord [m]

control_vane_length = 0.5 * duct_diameter  # one control vane length [m]
control_vane_chord = 0.30  # control vane chord [m]

propulsor_type = "traditional"  # options are: traditional, hybrid
n_blades = 6  # number of propeller blades [-]
rpm = 1000  # RPM of the propulsor
c_root = 0.2  # root chord of the propeller blade [m]
c_tip = 0.2  # tip chord of the propeller blade [m]
hub_diameter = 0.6  # spinner hub diameter [m]
propeller_sweep = 0  # propeller sweep [deg]
propeller_pitch = 0  # propeller pitch [deg]

""" Airfoil profiles used -> from the NACA 4-series"""
pylon_airfoil = "0012"
duct_airfoil = "0012"
control_vanes_airfoil = "0016"
support_airfoil = "0012"
prop_airfoil = "Hamilton568F"

d_exit = 1.05 * duct_diameter

""" Test condition"""
power_condition = "on"

n_pax = 68  # number of passengers in the aircraft [-]
w_pax = 80 + 22  # 80 kg weight of passenger and 22 kg of lugage [kg]
