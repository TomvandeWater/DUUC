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
from analysis_modules.plotting_functions import (plot_inflow_properties2, plot_weight_distribution, weight_comp_table,
                                                 interference_drag_range, cd0_drag_comparison, cd0_drag_empennage,
                                                 cl_cd_bucket,
                                                 cd_interference_drag_comparison, cl_comparison, xplot)
from analysis_modules.htail_sizing import slopes
from analysis_modules.BEM.BEM_sequence import bem_sequence
import matlab.engine
BEM_matlab_engine = matlab.engine.start_matlab()


""" Initialization of the analysis module """
alpha = 0

v_cruise = 128  # cruise speed [m/s] circa 250 knts
v_array = np.linspace(v_cruise - 1, v_cruise + 1, 3)

density_cruise = air_density_isa(flow_conditions.altitude)[0]
temperature_cruise = air_density_isa(flow_conditions.altitude)[2]

advance_duuc = advance_ratio(v_cruise, config.rpm, config.duct_diameter)
print(f"advance DUUC: {advance_duuc}")
advance_atr = advance_ratio(v_cruise, config.rpm, ref.blade_diameter)

""" RUN BEM model for both aircrafts"""
BEM_matlab_engine.cd(r'C:\Users\tomva\pythonProject\DUUC\analysis_modules\BEM')

# thrustdata = bem_sequence(config.n_blades, config.duct_diameter, v_cruise, 0.1, 1.2, density_cruise, temperature_cruise, 'HM568F')
# print(f"Thrust vector: {thrustdata}")

duuc_t_out, duuc_q_out, duuc_n_out, duuc_tc, duuc_cp, duuc_ct, duuc_va, duuc_vt = BEM_matlab_engine.BEM2(config.n_blades,
                                                                                       config.duct_diameter,
                                                                                       25, 0, v_cruise, 0, advance_duuc,
                                                                                       density_cruise,
                                                                                       temperature_cruise,
                                                                                       'HM568F', nargout=8)
bem_duuc = [duuc_t_out, duuc_q_out, duuc_n_out, duuc_tc, duuc_cp, duuc_ct, duuc_va, duuc_vt]

atr_t_out, atr_q_out, atr_n_out, atr_tc, atr_cp, atr_ct, atr_va, atr_vt = BEM_matlab_engine.BEM2(config.n_blades, ref.blade_diameter,
                                                                                 10, 0, v_cruise, 0, advance_atr,
                                                                                 density_cruise, temperature_cruise,
                                                                                 'HM568F', nargout=8)
bem_atr = [atr_t_out, atr_q_out, atr_n_out, atr_tc, atr_cp, atr_ct, atr_va, atr_vt]
BEM_matlab_engine.quit()

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
    power_condition = "on"

    """ Calculate mach number """
    mach_cal = v_inf / (speed_of_sound(flow_conditions.altitude))  # free stream mach number [-]
    mach_cal = 0.44
    """ Calculate Reynolds number"""
    density = air_density_isa(flow_conditions.altitude)

    reynolds_wing = reynolds(density, v_inf, ref.c_mac_w)
    #print(reynolds_wing)
    reynolds_duct = reynolds(density, v_inf, config.duct_chord)
    reynolds_htail = reynolds(density, v_inf, ref.c_root_h)

    advance = advance_ratio(v_inf, config.rpm, config.duct_diameter)

    """ setup the DUUC object from propulsive empennage class"""
    duuc: Aircraft = Aircraft(aircraft_type="DUUC", alpha=alpha, v_inf=v_inf, mach=mach_cal,
                              power=power_condition, bem_input=bem_duuc, delta_e=flow_conditions.delta_e,
                              delta_r=flow_conditions.delta_r, altitude=flow_conditions.altitude,
                              rpm=config.rpm,
                              n_blades=config.n_blades,
                              prop_diameter=0.95*config.duct_diameter,
                              prop_airfoil=config.prop_airfoil,
                              prop_c_root=config.c_root,
                              prop_c_tip=config.c_tip,
                              d_duct=config.duct_diameter,
                              c_duct=config.duct_chord,
                              duct_airfoil=config.duct_airfoil,
                              cant_angle=config.cant_angle,
                              b_pylon=config.pylon_length,
                              c_pylon=config.pylon_chord,
                              pylon_airfoil=config.pylon_airfoil,
                              l_nacelle=config.nacelle_length,
                              d_nacelle=config.nacelle_diameter,
                              b_support=config.support_length,
                              c_support=config.support_chord,
                              support_airfoil=config.support_airfoil,
                              l_cv=config.control_vane_length,
                              c_cv=config.control_vane_chord,
                              cv_airfoil=config.control_vanes_airfoil,
                              propulsor_type=config.propulsor_type,
                              d_hub=config.hub_diameter,
                              fus_d=ref.diameter_fuselage,
                              l_cab=ref.l_cab,
                              l_coc=ref.l_cockpit,
                              l_tail=ref.l_tail,
                              pax=config.n_pax,
                              wing_airfoil=ref.wing_airfoil,
                              wing_cr=ref.c_root_w,
                              wing_tr=ref.tr_w,
                              wing_span=ref.b_w,
                              wing_sweep=ref.phi_qc_w)

    """ data structure: --  v_inf|Duct|prop_f|prop_af|support|cv|wake|free_stream  """
    v_inflow = [v_inf, duuc.empennage.duct.inflow_velocity(), duuc.empennage.propeller.inflow_velocity(),
                duuc.empennage.support.inflow_velocity(), duuc.empennage.support.inflow_velocity(),
                duuc.empennage.elevator.inflow_velocity(), duuc.empennage.propeller.u1(), v_inf]

    a_inflow = [alpha, duuc.empennage.duct.inflow_angle(), duuc.empennage.propeller.inflow_angle(),
                duuc.empennage.support.inflow_angle(), duuc.empennage.support.inflow_angle(),
                duuc.empennage.elevator.inflow_angle(), flow_conditions.delta_e / 2, alpha]

    w_duuc = [duuc.fuselage.weight(), duuc.wing.weight(), duuc.empennage.duct.weight(),
              duuc.empennage.weight_ps()[0], duuc.empennage.weight_ps()[1],
              duuc.empennage.elevator.weight(), duuc.empennage.propeller.weight_engine(),
              duuc.empennage.nacelle.weight(), duuc.empennage.propeller.weight_fan()]

    station = [0, 1, 1.8, 2.2, 2.55, 4.5, 5.5, 6]
    # plot_inflow_properties2(v_inflow, a_inflow, station)
    # print(f"v inflow: {v_inflow}")
    # print(f"a inflow: {a_inflow}")

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
    atr: Aircraft = Aircraft(aircraft_type="conventional", alpha=alpha, v_inf=v_inf, mach=mach_cal,
                             power=power_condition, bem_input=bem_atr, delta_e=flow_conditions.delta_e,
                             delta_r=flow_conditions.delta_r, altitude=flow_conditions.altitude,
                             rpm=config.rpm,
                             n_blades=config.n_blades,
                             prop_diameter=0.95 * config.duct_diameter,
                             prop_airfoil=config.prop_airfoil,
                             prop_c_root=config.c_root,
                             prop_c_tip=config.c_tip,
                             d_duct=config.duct_diameter,
                             c_duct=config.duct_chord,
                             duct_airfoil=config.duct_airfoil,
                             cant_angle=config.cant_angle,
                             b_pylon=config.pylon_length,
                             c_pylon=config.pylon_chord,
                             pylon_airfoil=config.pylon_airfoil,
                             l_nacelle=config.nacelle_length,
                             d_nacelle=config.nacelle_diameter,
                             b_support=config.support_length,
                             c_support=config.support_chord,
                             support_airfoil=config.support_airfoil,
                             l_cv=config.control_vane_length,
                             c_cv=config.control_vane_chord,
                             cv_airfoil=config.control_vanes_airfoil,
                             propulsor_type=config.propulsor_type,
                             d_hub=config.hub_diameter,
                             fus_d=ref.diameter_fuselage,
                             l_cab=ref.l_cab,
                             l_coc=ref.l_cockpit,
                             l_tail=ref.l_tail,
                             pax=config.n_pax,
                             wing_airfoil=ref.wing_airfoil,
                             wing_cr=ref.c_root_w,
                             wing_tr=ref.tr_w,
                             wing_span=ref.b_w,
                             wing_sweep=ref.phi_qc_w)

    w_atr = [atr.fuselage.weight(), atr.wing.weight(), atr.empennage.propeller.weight_engine(),
             atr.empennage.ht_tail.weight(), atr.empennage.vt_tail.weight(), atr.empennage.weight_cv(),
             atr.empennage.weight_nac() / 2, atr.empennage.propeller.weight_engine() / 2,
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
    cl_nac = atr.empennage.nacelle.cl() * 2
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

#print(f"S_wet atr: {atr.empennage.area_wet()}")
#print(f"S_wet duuc: {duuc.empennage.area_wet()}")
# print(results_matrix)

""" Plot section """
weight_comp_table(w_atr, "conventional")
weight_comp_table(w_duuc, "DUUC")
# plot_weight_distribution(w_duuc, w_atr)

# cd0_drag_comparison(atr.cd0_vector(), duuc.cd0_vector())
# print(f"cd0 atr: {sum(atr.cd0_empennage())}, cd0 duuc: {sum(duuc.cd0_empennage())}")


cd0_drag_empennage(atr.cd0_empennage(), duuc.cd0_empennage())

""" Polar of the ATR with reference"""
# cdi_cal = 1 / (np.pi * atr.wing.aspect_ratio() * atr.wing.oswald())
# cl_cd_bucket(sum(atr.cd0_vector()), cdi_cal, alpha)

""" Plot to compare interference drag per component"""
# cd_interference_drag_comparison(atr.empennage.cd_interference_vector(), duuc.empennage.cd_interference_vector())
# interference_drag_range(cdi_results)

""" Plot to compare CL de-normalized per component"""
# cl_comparison(atr.cl_vector(), duuc.cl_vector(), alpha)
duuc.x_cog()
atr.x_cog()

""" Horizontal tail sizing """
# xplot(-0.4, 0.189, 0.20, (atr.x_cog()[0] - atr.x_cog()[1] - 0.25), 5, ref.c_mac_w, 0.25 * ref.c_mac_w, "ATR", bem_atr[0])
# xplot(duuc_a1, duuc_b1, duuc_a2, cog_duuc - x_lemac_duuc, 5, ref.c_mac_w, 0.25 * ref.c_mac_w, "DUUC", bem_duuc[0])
