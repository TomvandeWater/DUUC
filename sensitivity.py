import matplotlib.pyplot as plt
import config
import data.atr_reference as ref
import numpy as np
from analysis_modules.aerodynamic import *
from analysis_modules.factors import *
from aircraft.aircraft_assembly import Aircraft

parameters = {
            'duct_diameter': config.duct_diameter, 'duct_chord': config.duct_chord, 'duct_profile': config.duct_airfoil,

            'pylon_chord': config.pylon_chord, 'pylon_length': config.pylon_length, 'cant_angle': config.cant_angle,
            'pylon_profile': config.pylon_airfoil,

            'support_chord': config.support_chord, 'support_length': config.support_length,
            'support_profile': config.support_airfoil,

            'num_blades': config.n_blades, 'propeller_diameter': config.duct_diameter - 0.22,
            'hub_diameter': config.hub_diameter, 'propeller_airfoil': config.prop_airfoil,
            'propeller_c_root': config.c_root, 'propeller_c_tip': config.c_tip,

            'BEM1': 4142, 'BEM2': 2648, 'BEM4': -1.44, 'BEM5': 0.889,
            'BEM6': 0.329, 'BEM7': 5, 'BEM8': 10, 'BEM3': 1820,

            'altitude': 7000,  'velocity': 128, 'alpha': 0.0, 'delta_e': 0, 'delta_r': 0,
            'power_condition': 'on', 'propulsion_type': 'conventional', 'RPM': config.rpm,
            "aircraft_n_pax": config.n_pax, "static_margin": 5, "beta": 0,

            'wing_span': ref.b_w, 'wing_phi_qc': ref.phi_qc_w, 'wing_airfoil': ref.wing_airfoil, 'wing_tr': ref.tr_w,
            'wing_c_root': np.round(ref.c_root_w, 3),

            "fuselage_length": np.round((ref.l_tail+ref.l_cockpit+ref.l_cabin), 3),
            "fuselage_diameter": ref.diameter_fuselage, "fuselage_co_l": ref.l_cockpit, "fuselage_ca_l": ref.l_cabin,
            "fuselage_ta_l": ref.l_tail,

            "nacelle_length": config.nacelle_length, "nacelle_diameter": config.nacelle_diameter,

            "hcv_span": config.control_vane_length, "hcv_chord": config.control_vane_chord,
            "vcv_span": config.control_vane_length, "vcv_chord": config.control_vane_chord,
            "cv_airfoil": config.control_vanes_airfoil,

            "l_v": 9.13, "x_prop": 0.3, "y_engine": ref.y_engine, "x_support": 0.5 * config.duct_chord,
            "x_control_vanes": 0.95 * config.duct_chord, "x_wing": 11.5, "x_pylon": 0.5 * config.duct_chord,
            "a_i_wing": 0, "a_i_duct": 0, "x_PE": 30, "y_PE": 0, # -> y parameter is unused right now as
            "z_PE": ref.diameter_fuselage,
        }

duct_diameter = parameters.get('duct_diameter')
duct_chord = parameters.get('duct_chord')
duct_profile = parameters.get('duct_profile')

pylon_chord = parameters.get('pylon_chord')
pylon_length = parameters.get('pylon_length')

support_chord = parameters.get('support_chord')
support_length = parameters.get('support_length')

num_blades = parameters.get('num_blades')
altitude = parameters.get('altitude')
velocity = parameters.get('velocity')
alpha = parameters.get('alpha')
mach = velocity / speed_of_sound(altitude)
delta_e = parameters.get('delta_e')
delta_r = parameters.get('delta_r')
beta = parameters.get('beta')
power_condition = parameters.get('power_condition')
propulsion_type = parameters.get('propulsion_type')

cant_angle = parameters.get('cant_angle')
pylon_profile = parameters.get('pylon_profile')
support_profile = parameters.get('support_profile')
BEM1 = parameters.get('BEM1')
BEM2 = parameters.get('BEM2')
BEM3 = parameters.get('BEM3')
BEM4 = parameters.get('BEM4')
BEM5 = parameters.get('BEM5')
BEM6 = parameters.get('BEM6')
BEM7 = parameters.get('BEM7')
BEM8 = parameters.get('BEM8')
BEM9 = parameters.get('BEM9')

bem_input = [BEM1, BEM2, BEM3, BEM4, BEM5, BEM6, BEM7, BEM8, BEM9]
RPM = parameters.get('RPM')
propeller_diameter = parameters.get('propeller_diameter')
hub_diameter = parameters.get('hub_diameter')
propeller_airfoil = parameters.get('propeller_airfoil')
propeller_c_root = parameters.get('propeller_c_root')
propeller_c_tip = parameters.get('propeller_c_tip')

wing_span = parameters.get('wing_span')
wing_phi_qc = parameters.get('wing_phi_qc')
wing_airfoil = parameters.get('wing_airfoil')
wing_tr = parameters.get('wing_tr')
wing_c_root = parameters.get('wing_c_root')
x_wing = parameters.get('x_wing')

fuselage_diameter = parameters.get('fuselage_diameter')
fuselage_co_l = parameters.get('fuselage_co_l')
fuselage_ca_l = parameters.get('fuselage_ca_l')
fuselage_ta_l = parameters.get('fuselage_ta_l')
fuselage_length = fuselage_ca_l + fuselage_ta_l + fuselage_co_l

aircraft_n_pax = parameters.get('aircraft_n_pax')
x_PE = parameters.get('x_PE')
z_PE = parameters.get('z_PE')
y_engine = parameters.get('y_engine')

nacelle_length = parameters.get('nacelle_length')
nacelle_diameter = parameters.get('nacelle_diameter')

hcv_span = parameters.get('hcv_span')
hcv_chord = parameters.get('hcv_chord')

cv_airfoil = parameters.get('cv_airfoil')
x_prop = parameters.get('x_prop') * parameters.get("duct_chord")
x_pylon = parameters.get('x_pylon')
x_support = parameters.get('x_support')
x_cv = parameters.get('x_cv')
a_install_wing = parameters.get('a_i_wing')
a_install_duct = parameters.get('a_i_duct')
stm = parameters.get('static_margin')

conditions = [velocity, alpha, altitude, mach, a_install_wing, a_install_duct, beta]
reference = [ref.s_w, ref.c_mac_w]
geometry_duct = [duct_diameter, duct_chord, duct_profile]
geometry_pylon = [pylon_length, pylon_chord, pylon_profile, cant_angle]
geometry_support = [support_length, support_chord, support_profile, cant_angle]
geometry_control = [hcv_span, hcv_chord, cv_airfoil]
geometry_propeller = [num_blades, propeller_diameter, hub_diameter, propeller_airfoil, 0, 0,
                      propeller_c_root, propeller_c_tip]
geometry_nacelle = [nacelle_length, nacelle_diameter]

comp_pe = [x_pylon, x_support, x_prop, x_cv]

duuc: Aircraft = Aircraft(aircraft_type="DUUC", conditions=conditions,
                          power=power_condition, bem_input=bem_input, delta_e=delta_e,
                          delta_r=delta_r, reference=reference,
                          wing_span=wing_span, wing_sweep=wing_phi_qc, wing_airfoil=wing_airfoil,
                          wing_tr=wing_tr, wing_cr=wing_c_root, l_coc=fuselage_co_l, l_cab=fuselage_ca_l,
                          l_tail=fuselage_ta_l, fus_d=fuselage_diameter, pax=aircraft_n_pax,
                          geometry_duct=geometry_duct, geometry_pylon=geometry_pylon,
                          geometry_control=geometry_control, geometry_support=geometry_support,
                          geometry_nacelle=geometry_nacelle, geometry_propeller=geometry_propeller,
                          geometry_ht=[], geometry_vt=[], x_duct=x_PE, x_wing=x_wing, comp_pe=comp_pe,
                          propulsor_type=propulsion_type, rpm=RPM, z_duct=z_PE)




