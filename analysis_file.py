import flow_conditions
import config
from analysis_modules.general_functions import shp_conversion
from analysis_modules.aerodynamic import reynolds, propulsive_efficiency, speed_of_sound
from analysis_modules.ISA import air_density_isa
from analysis_modules.factors import advance_ratio
import matplotlib.pyplot as plt
from aircraft.aircraft_assembly import Aircraft
import numpy as np
import data.atr_reference as ref
from analysis_modules.plotting_functions import plot_inflow_properties2, plot_weight_distribution, weight_comp_table, interference_drag_range


""" Initialization of the analysis module """
alpha = 0

v_cruise = 128  # cruise speed [m/s] circa 250 knts
v_array = np.linspace(v_cruise - 50, v_cruise + 20, 71)

""" setup matrices for result printing """
cd_results = []
cdi_results = []
cd0_results = []
cd_trad = []
cl_trad = []
cl_results = []
eta_prop_results = []
results = []
eta_prop_debug = []

for i in range(len(v_array)):
    v_inf = float(v_array[i])
    va_inlet = v_inf * (np.pi * (config.duct_diameter / 2))
    c = 9.81  # conversion factor from Newton to kg

    """ Calculate mach number """
    mach_cal = v_inf / (speed_of_sound(flow_conditions.altitude))  # free stream mach number [-]

    """ Calculate Reynolds number"""
    reynolds_wing = reynolds(air_density_isa(flow_conditions.altitude), v_inf, ref.c_root_w)

    advance = advance_ratio(v_inf, config.rpm, config.duct_diameter)

    """ setup the DUUC object from propulsive empennage class"""
    duuc: Aircraft = Aircraft(aircraft_type="DUUC", alpha=alpha, reynolds=reynolds_wing, v_inf=v_inf,
                              mach=mach_cal)

    """ data structure: --  v_inf|Duct|prop_f|prop_af|support|cv|wake|free_stream  """
    v_inflow = [v_inf, duuc.empennage.duct.inflow_velocity(), duuc.empennage.propeller.inflow_velocity(),
                duuc.empennage.support.inflow_velocity(), duuc.empennage.support.inflow_velocity(),
                duuc.empennage.elevator.inflow_velocity(), duuc.empennage.propeller.u1(), v_inf]

    a_inflow = [alpha, duuc.empennage.duct.inflow_angle(), duuc.empennage.propeller.inflow_angle(),
                duuc.empennage.support.inflow_angle(), duuc.empennage.support.inflow_angle(),
                duuc.empennage.elevator.inflow_angle(), flow_conditions.delta_e / 2, alpha]

    w_duuc = [duuc.fuselage.weight(), duuc.wing.weight(), duuc.empennage.duct.weight(),
              duuc.empennage.pylon.weight(), duuc.empennage.support.weight(),
              duuc.empennage.elevator.weight(), duuc.empennage.propeller.weight_engine(),
              duuc.empennage.nacelle.weight(), duuc.empennage.propeller.weight_fan()]

    station = [0, 1, 1.8, 2.2, 2.55, 4.5, 5.5, 6]
    # plot_inflow_properties2(v_inflow, a_inflow, station)

    """ cd properties """
    cd_total = duuc.empennage.cd_prime()
    cd_duct = duuc.empennage.duct.cd_prime()
    cd_pylon = duuc.empennage.pylon.cd_prime()
    cd_support = duuc.empennage.support.cd_prime()
    cd_control = 2 * duuc.empennage.elevator.cd_prime() + 2 * duuc.empennage.rudder.cd_prime()
    cd_nacelle = duuc.empennage.nacelle.cd_prime()
    cd_prop = duuc.empennage.propeller.cd_prime()

    cdi_duct = duuc.empennage.duct.cdi()
    cdi_pylon = duuc.empennage.pylon.cd_interference()
    cdi_support = duuc.empennage.support.cd_interference()
    cdi_control = 2 * duuc.empennage.elevator.cdi() + 2 * duuc.empennage.rudder.cdi()
    cdinter_control = 2 * duuc.empennage.elevator.cd_interference() + 2 * duuc.empennage.rudder.cd_interference()
    cdi_nacelle = 0
    cdi_prop = duuc.empennage.propeller.cd_interference()

    cd0_duct = duuc.empennage.duct.cd0()
    cd0_pylon = duuc.empennage.pylon.cd()
    cd0_support = duuc.empennage.support.cd()
    cd0_control = 2 * duuc.empennage.elevator.cd0() + 2 * duuc.empennage.rudder.cd0()
    cd0_nacelle = duuc.empennage.nacelle.cd0()
    cd0_prop = duuc.empennage.propeller.cd0()

    """ cl properties """
    cl_total = duuc.empennage.cl_prime()
    cl_duct = duuc.empennage.duct.cl_prime()
    cl_pylon = duuc.empennage.pylon.cl_prime()
    cl_support = duuc.empennage.support.cl_prime()
    cl_control = 2 * duuc.empennage.elevator.cl_prime() + 2 * duuc.empennage.rudder.cl_prime()
    cl_nacelle = duuc.empennage.nacelle.cl_prime()
    cl_prop = duuc.empennage.propeller.cl_prime()

    """ data appending from instance """
    cd_results.append([advance, cd_total, cd_duct, cd_pylon, cd_support, cd_control, cd_nacelle, cd_prop])
    cdi_results.append([v_inf, cdinter_control, cdi_duct, cdi_pylon, cdi_support, cdi_control, cdi_nacelle, cdi_prop])
    cd0_results.append([v_inf, "cd0_total", cd0_duct, cd0_pylon, cd0_support, cd0_control, cd0_nacelle, cd0_prop])
    cl_results.append([v_inf, cl_total, cl_duct, cl_pylon, cl_support, cl_control, cl_nacelle, cl_prop])

    """ define reference aircraft class """
    atr: Aircraft = Aircraft(aircraft_type="conventional", alpha=alpha, reynolds=reynolds_wing, v_inf=v_inf,
                             mach=mach_cal)

    w_atr = [atr.fuselage.weight(), atr.wing.weight(), atr.empennage.propeller.weight_engine(),
             atr.empennage.ht_tail.weight(), atr.empennage.vt_tail.weight(), atr.empennage.weight_cv(),
             atr.empennage.weight_nac() / 2, atr.empennage.propeller.weight_engine(),
             atr.empennage.propeller.weight_fan()]

    """ cd properties """
    cd_total_trad = atr.empennage.cd_prime()
    cd_vt = atr.empennage.vt_tail.cd_prime()
    cd_ht = atr.empennage.ht_tail.cd_prime()
    cd_nac = atr.empennage.nacelle.cd_prime() * 2
    cd_prop_trad = atr.empennage.propeller.cd_prime() * 2

    """ cl properties """
    cl_total_trad = atr.empennage.cl_prime()
    cl_ht = atr.empennage.ht_tail.cl_prime()
    cl_nac = atr.empennage.nacelle.cl_prime() * 2
    cl_prop_trad = atr.empennage.propeller.cl_prime() * 2

    """ data appending from instance """
    cd_trad.append([v_inf, cd_total_trad, cd_vt, cd_ht, cd_nac, cd_prop_trad])
    cl_trad.append([v_inf, cl_total_trad, 0, cl_ht, cl_nac, cl_prop_trad])

    """ determine propulsive efficiency """
    fx_duuc = duuc.thrust() - duuc.empennage.drag()
    power_in_duuc = shp_conversion(ref.P_max_cruise, "SI")

    fx_ref = atr.thrust() - atr.empennage.drag()
    power_in_ref = shp_conversion(ref.P_max_cruise, "SI")

    """ values for recovery"""
    k_srv = 0.97
    k_tip = 0.95

    """ calculate propulsive efficiency of both empennages """
    eta_prop_duuc = propulsive_efficiency(fx_duuc, v_inf, power_in_duuc, k_srv, k_tip)
    eta_prop_ref = propulsive_efficiency(fx_ref, v_inf, power_in_ref, k1=1, k2=1)

    eta_prop_results.append([advance, eta_prop_duuc, eta_prop_ref])
    eta_prop_debug.append([v_inf, power_in_duuc, power_in_ref, duuc.thrust(), atr.thrust(), duuc.empennage.drag(),
                           atr.empennage.drag()])
    results.append([v_inf, duuc.thrust(), duuc.empennage.drag(), fx_duuc, "|", atr.thrust(), atr.empennage.drag(),
                    fx_ref, '|', power_in_duuc, power_in_ref])

""" convert to numpy array for manipulation """
cd_results_conv_matrix = np.array(cd_trad)
cl_results_conv_matrix = np.array(cl_trad)
cl_results_matrix = np.array(cl_results)
results_matrix = np.array(results)
eta_debug = np.array(eta_prop_debug)

# print(results_matrix)

""" Plot section """
weight_comp_table(w_atr, "conventional")
weight_comp_table(w_duuc, "DUUC")
plot_weight_distribution(w_duuc, w_atr)
interference_drag_range(cdi_results)
