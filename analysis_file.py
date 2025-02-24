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
from analysis_modules.plotting_functions import plot_inflow_velocity, plot_inflow_angle


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

    """ Calculate mach number """
    mach_cal = v_inf / (speed_of_sound(flow_conditions.altitude))  # free stream mach number [-]

    """ Calculate Reynolds number"""
    reynolds_wing = reynolds(air_density_isa(flow_conditions.altitude),
                             v_inf, ref.c_root_w)
    # print(f"Re duct: {duct_reynolds}")

    advance = advance_ratio(v_inf, config.rpm, config.duct_diameter)

    """ setup the DUUC object from propulsive empennage class"""
    duuc: Aircraft = Aircraft(aircraft_type="DUUC", alpha=alpha, reynolds=reynolds_wing, v_inf=v_inf,
                              mach=mach_cal)

    v_inflow = [v_inf, duuc.empennage.duct.inflow_velocity(), duuc.empennage.propeller.inflow_velocity(),
                duuc.empennage.support.inflow_velocity(),
                duuc.empennage.elevator.inflow_velocity(), duuc.empennage.propeller.u1(), 100, 100]
    a_inflow = [alpha, duuc.empennage.duct.inflow_angle(), duuc.empennage.propeller.inflow_angle(),
                duuc.empennage.support.inflow_angle(),
                duuc.empennage.elevator.inflow_angle(), 0, 0, 0]

    station = [0, 1, 1.8, 2.2, 2.55, 4.5, 5.5, 6]

    # plot_inflow_velocity(v_inflow, a_inflow, station)
    # plot_inflow_angle(a_inflow, station)

    """ cd properties """
    cd_total = duuc.empennage.cd_prime()
    print(f"DUUC cdprime: {cd_total}")
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
    eta_prop_debug.append([v_inf, power_in_duuc, power_in_ref, duuc.thrust(), atr.thrust(), duuc.empennage.drag(), atr.empennage.drag()])
    results.append([v_inf, duuc.thrust(), duuc.empennage.drag(), fx_duuc, "|", atr.thrust(), atr.empennage.drag(), fx_ref, '|', power_in_duuc, power_in_ref])

""" convert to numpy array for manipulation """
cd_results_matrix = np.array(cd_results)
cdi_results_matrix = np.array(cdi_results)
cd0_results_matrix = np.array(cd0_results)
cd_results_conv_matrix = np.array(cd_trad)
cl_results_conv_matrix = np.array(cl_trad)
cl_results_matrix = np.array(cl_results)
eta_results_matrix = np.array(eta_prop_results)
results_matrix = np.array(results)
eta_debug = np.array(eta_prop_debug)

# print(results_matrix)

""" Plot section """
"""---------           plot about drag components in the DUUC                           ---------"""
plt.figure('Drag polar DUUC')
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 1], label=r'Total drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 2], label=r'Duct drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 3], label=r'Pylon drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 4], label=r'Support drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 5], label=r'Control drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 6], label=r'Nacelle drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 7], label=r'Propeller drag coefficient')
plt.xlabel(r'$v_{inf}$')
plt.ylabel(r'$C_{d}$')
plt.title(r'$C_{d}$ polar of duct')
plt.legend()
plt.grid(True)

"""---------           plot about drag components in conventional empennage             ---------"""
plt.figure('Drag polar ATR 72')
plt.plot(cd_results_conv_matrix[:, 0], cd_results_conv_matrix[:, 1], label=r'Total drag coefficient')
plt.plot(cd_results_conv_matrix[:, 0], cd_results_conv_matrix[:, 2], label=r'VT drag coefficient')
plt.plot(cd_results_conv_matrix[:, 0], cd_results_conv_matrix[:, 3], label=r'HT drag coefficient')
plt.plot(cd_results_conv_matrix[:, 0], cd_results_conv_matrix[:, 4], label=r'Nacelle drag coefficient')
plt.plot(cd_results_conv_matrix[:, 0], cd_results_conv_matrix[:, 5], label=r'Propeller drag coefficient')
plt.xlabel(r'$v_{inf}$')
plt.ylabel(r'$C_{d}$')
plt.title(r'$C_{d}$ polar of reference aircraft')
plt.legend()
plt.grid(True)

"""---------           plot about lift components in DUUC                               ---------"""
plt.figure('Lift polar duct')
plt.plot(cl_results_matrix[:, 0], cl_results_matrix[:, 1], label=r'Total lift coefficient')
plt.plot(cl_results_matrix[:, 0], cl_results_matrix[:, 2], label=r'Duct lift coefficient')
plt.plot(cl_results_matrix[:, 0], cl_results_matrix[:, 3], label=r'Pylon lift coefficient')
plt.plot(cl_results_matrix[:, 0], cl_results_matrix[:, 4], label=r'Support lift coefficient')
plt.plot(cl_results_matrix[:, 0], cl_results_matrix[:, 5], label=r'Control lift coefficient')
plt.plot(cl_results_matrix[:, 0], cl_results_matrix[:, 6], label=r'Nacelle lift coefficient')
plt.plot(cl_results_matrix[:, 0], cl_results_matrix[:, 7], label=r'Propeller lift coefficient')
plt.xlabel(r'$v_{inf}$')
plt.ylabel(r'$C_{l}$')
plt.title(r'$C_{l}$ polar of duct')
plt.legend()
plt.grid(True)

"""---------           plot about propulsive efficiency             ---------"""
plt.figure('Propulsive efficiency')
plt.plot(eta_results_matrix[:, 0], eta_results_matrix[:, 1], label=r'DUUC')
plt.plot(eta_results_matrix[:, 0], eta_results_matrix[:, 2], label=r'ATR72-600')
plt.xlabel(r'J [-]')
plt.ylabel(r'$\eta_{prop}$ [%]')
plt.title(r'Propulsive efficiency')
plt.legend()
plt.grid(True)

"""---------            plot cl-cd bucket                           ----------"""
plt.figure('Cl-cd bucket')
plt.plot(cd_results_matrix[:, 1] + cd_results_matrix[:, 3], cl_results_matrix[:, 1] + cl_results_matrix[:, 3],label=r'DUUC')
plt.plot(cd_results_conv_matrix[:, 1], cl_results_conv_matrix[:, 1], label=r'ATR72-600')
plt.xlabel(r'$C_{d}$ [-]')
plt.ylabel(r'$C_{l}$ [-]')
plt.title(r'Comparison of DUUC and ATR72')
plt.legend()
plt.grid(True)
"""---------           debug plot about propulsive efficiency             ---------"""
""" stays constant -> checked
plt.figure('Power in')
plt.plot(eta_debug[:, 0], eta_debug[:, 1], label=r'Power in DUUC')
plt.plot(eta_debug[:, 0], eta_debug[:, 2], label=r'Power in ATR72-600')
plt.xlabel(r'$v_{inf} [m/s]$')
plt.ylabel(r'$Power [W]$')
plt.title(r'Power in')
plt.legend()
plt.grid(True)"""

plt.figure('Force comparison')
plt.plot(eta_debug[:, 0], eta_debug[:, 3], label=r'Thrust DUUC')
plt.plot(eta_debug[:, 0], eta_debug[:, 4], label=r'Thrust ATR72-600')
plt.plot(eta_debug[:, 0], eta_debug[:, 5], label=r'Drag DUUC')
plt.plot(eta_debug[:, 0], eta_debug[:, 6], label=r'Drag ATR72-600')
plt.xlabel(r'$v_{inf}$ [m/s]')
plt.ylabel(r'$Force$ [N]')
plt.title(r'Forces')
plt.legend()
plt.grid(True)

""" Drag figures"""
"""
plt.figure("Duct drag coefficients")
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 2], label=r'Total drag coefficient')
plt.plot(cd_results_matrix[:, 0], cdi_results_matrix[:, 2], label=r'Induced drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd0_results_matrix[:, 2], label=r'Zero lift drag coefficient')
plt.xlabel(r'$v_{inf} [m/s]$')
plt.ylabel(r'$C_{d} [-]$')
plt.title(r'Duct drag coefficients')
plt.legend()

plt.figure("Pylon drag coefficients")
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 3], label=r'Total drag coefficient')

plt.plot(cd_results_matrix[:, 0], cd0_results_matrix[:, 3], label=r'Zero lift drag coefficient')
plt.xlabel(r'$v_{inf} [m/s]$')
plt.ylabel(r'$C_{d} [-]$')
plt.title(r'Pylon drag coefficients')
plt.legend()


plt.figure("Support drag coefficients")
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 4], label=r'Total drag coefficient')

plt.plot(cd_results_matrix[:, 0], cd0_results_matrix[:, 4], label=r'Zero lift drag coefficient')
plt.xlabel(r'$v_{inf} [m/s]$')
plt.ylabel(r'$C_{d} [-]$')
plt.title(r'Support drag coefficients')
plt.legend()

plt.figure("Control Vanes drag coefficients")
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 5], label=r'Total drag coefficient')
plt.plot(cd_results_matrix[:, 0], cdi_results_matrix[:, 5], label=r'Induced drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd0_results_matrix[:, 5], label=r'Zero lift drag coefficient')
plt.xlabel(r'$v_{inf} [m/s]$')
plt.ylabel(r'$C_{d} [-]$')
plt.title(r'Control Vanes drag coefficients')
plt.legend()


plt.figure("Nacelle drag coefficients")
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 6], label=r'Total drag coefficient')
plt.plot(cd_results_matrix[:, 0], cdi_results_matrix[:, 6], label=r'Induced drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd0_results_matrix[:, 6], label=r'Zero lift drag coefficient')
plt.xlabel(r'$v_{inf} [m/s]$')
plt.ylabel(r'$C_{d} [-]$')
plt.title(r'Nacelle drag coefficients')
plt.legend()

plt.figure("Propeller drag coefficients")
plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, 7], label=r'Total drag coefficient')
plt.plot(cd_results_matrix[:, 0], cdi_results_matrix[:, 7], label=r'Interference drag coefficient')
plt.plot(cd_results_matrix[:, 0], cd0_results_matrix[:, 7], label=r'Zero lift drag coefficient')
plt.xlabel(r'$v_{inf} [m/s]$')
plt.ylabel(r'$C_{d} [-]$')
plt.title(r'Propeller drag coefficients')
plt.legend()"""


plt.figure("Interference drag coefficients")
plt.plot(cd_results_matrix[:, 0], cdi_results_matrix[:, 7], label=r'Propeller')
plt.plot(cd_results_matrix[:, 0], cdi_results_matrix[:, 1], label=r'Control Vanes')
plt.plot(cd_results_matrix[:, 0], cdi_results_matrix[:, 4], label=r'Support')
plt.plot(cd_results_matrix[:, 0], cdi_results_matrix[:, 3], label=r'Pylon')
plt.xlabel(r'J [-]')
plt.ylabel(r'$C_{d}$ [-]')
plt.title(r'Interference drag coefficients')
plt.grid(True)
plt.legend()

plt.show()
