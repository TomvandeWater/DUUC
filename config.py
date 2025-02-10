""" Propulsive empennage set up"""
duct_diameter = 3.6  # duct diameter in [m]
duct_chord = 1.5  # duct chord in [m]

pylon_length = 1  # pylon length in [m]
pylon_chord = 0.5  # pylon chord length [m]
cant_angle = 30  # cant angle of the pylon wrt to the fuselage [deg]

nacelle_length = 1.2  # nacelle length in [m]
nacelle_diameter = 0.30  # nacelle diameter in [m]

support_length = duct_diameter  # support length in [m]
support_chord = 0.5  # support chord [m]

control_vane_length = 0.5 * duct_diameter  # one control vane length [m]
control_vane_chord = 0.20  # control vane chord [m]

""" Airfoil profiles used -> from the NACA 4-series"""
pylon_airfoil = "0012"
duct_airfoil = "0012"
control_vanes_airfoil = "0016"
support_airfoil = "0012"
