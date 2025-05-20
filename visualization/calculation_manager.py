from aircraft.aircraft_assembly import Aircraft
from analysis_modules.aerodynamic import speed_of_sound
from analysis_modules.vtail_sizing import *
from analysis_modules.htail_sizing import *


def calculation_manager(parameters):
    # -----                 PARAMETER UNPACKING                 ----- #
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
    cv_mode = parameters.get('control_vane_mode')

    conditions = [velocity, alpha, altitude, mach, a_install_wing, a_install_duct, beta]
    reference = [ref.s_w, ref.c_mac_w]
    geometry_duct = [duct_diameter, duct_chord, duct_profile]
    geometry_pylon = [pylon_length, pylon_chord, pylon_profile, cant_angle]
    geometry_support = [support_length, support_chord, support_profile, cant_angle]
    geometry_control = [hcv_span, hcv_chord, cv_airfoil]
    geometry_propeller = [num_blades, propeller_diameter, hub_diameter, propeller_airfoil, 0, 0,
                          propeller_c_root, propeller_c_tip]
    geometry_nacelle = [nacelle_length, nacelle_diameter]
    geometry_ht = [ref.b_h, ref.c_root_h, ref.airfoil_ht, ref.phi_qc_h, ref.tr_h, ref.c_root_h]
    geometry_vt = [ref.b_v, ref.c_root_v, ref.airfoil_vt, ref.phi_qc_v, ref.tr_v, ref.c_tip_v]

    comp_pe = [x_pylon, x_support, x_prop, x_cv]

    # -----                     CREATE DUUC INSTANCE                    ------ #
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
                              propulsor_type=propulsion_type, rpm=RPM, z_duct=z_PE, cv_mode=cv_mode)

    alpha_vector = np.linspace(0, 15, 31)
    pe_cl_vector = []
    pe_cd_vector = []
    pe_cm_vector = []
    pylon_vector = []
    support_vector = []
    duct_vector = []
    control_vector = []
    nacelle_vector = []
    for r in range(len(alpha_vector)):
        conditions_vector = [velocity, alpha_vector[r], altitude, mach, a_install_wing, a_install_duct, beta]
        duuc2: Aircraft = Aircraft(aircraft_type="DUUC", conditions=conditions_vector,
                                   power=power_condition, bem_input=bem_input, delta_e=delta_e,
                                   delta_r=delta_r, reference=reference,
                                   wing_span=wing_span, wing_sweep=wing_phi_qc, wing_airfoil=wing_airfoil,
                                   wing_tr=wing_tr, wing_cr=wing_c_root, l_coc=fuselage_co_l, l_cab=fuselage_ca_l,
                                   l_tail=fuselage_ta_l, fus_d=fuselage_diameter, pax=aircraft_n_pax,
                                   geometry_duct=geometry_duct, geometry_pylon=geometry_pylon,
                                   geometry_control=geometry_control, geometry_support=geometry_support,
                                   geometry_nacelle=geometry_nacelle, geometry_propeller=geometry_propeller,
                                   geometry_ht=[], geometry_vt=[], x_duct=x_PE, x_wing=x_wing, comp_pe=comp_pe,
                                   propulsor_type=propulsion_type, rpm=RPM, z_duct=z_PE, cv_mode=cv_mode)
        pe_cl_vector.append(duuc2.empennage.cl_sum())
        pe_cd_vector.append(duuc2.empennage.cd_sum())
        pe_cm_vector.append(duuc2.empennage.cm_emp())
        duct_vector.append([2 * duuc2.empennage.duct.cl()[0], 2 * duuc2.empennage.duct.cd()[0], 2 * duuc2.empennage.duct.cm()[0],
                            2 * duuc2.empennage.duct.cl()[1], 2 * duuc2.empennage.duct.cd()[1], 2 * duuc2.empennage.duct.cm()[1]])
        pylon_vector.append([2 * duuc2.empennage.pylon.cl()[0], 2 * duuc2.empennage.pylon.cd()[0], 2 * duuc2.empennage.pylon.cm()[0],
                             2 * duuc2.empennage.pylon.cl()[1], 2 * duuc2.empennage.pylon.cd()[1], 2 * duuc2.empennage.pylon.cm()[1]])
        support_vector.append([2 * duuc2.empennage.support.cl()[0], 2 * duuc2.empennage.support.cd()[0],
                               2 * duuc2.empennage.support.cm()[0], 2 * duuc2.empennage.support.cl()[1],
                               2 * duuc2.empennage.support.cd()[1], 2 * duuc2.empennage.support.cm()[1]])
        control_vector.append([4 * duuc2.empennage.elevator.cl()[0],  4 * duuc2.empennage.rudder.cd()[0] +
                               4 * duuc2.empennage.elevator.cd()[0], 0, 4 * duuc2.empennage.elevator.cl()[1],
                               4 * duuc2.empennage.rudder.cd()[1] + 4 * duuc2.empennage.elevator.cd()[1], 0])
        nacelle_vector.append([2 * duuc2.empennage.nacelle.cl()[0], 2 * duuc2.empennage.nacelle.cd()[0],
                               2 * duuc2.empennage.nacelle.cm()[0], 2 * duuc2.empennage.nacelle.cl()[1],
                               2 * duuc2.empennage.nacelle.cd()[1], 2 * duuc2.empennage.nacelle.cm()[1]])

    pylon_vect = np.array(pylon_vector)
    duct_vect = np.array(duct_vector)
    support_vect = np.array(support_vector)
    control_vect = np.array(control_vector)
    nacelle_vect = np.array(nacelle_vector)

    geometry_propeller_conv = [num_blades, propeller_diameter, hub_diameter, propeller_airfoil, 0, 0,
                               propeller_c_root, propeller_c_tip]
    geometry_nacelle_conv = [nacelle_length, nacelle_diameter]

    x_vector = np.linspace(0, fuselage_length, 100)
    x_cg_vector = []
    for length in range(len(x_vector)):
        duuc3: Aircraft = Aircraft(aircraft_type="DUUC", conditions=conditions,
                                   power=power_condition, bem_input=bem_input, delta_e=delta_e,
                                   delta_r=delta_r, reference=reference,
                                   wing_span=wing_span, wing_sweep=wing_phi_qc, wing_airfoil=wing_airfoil,
                                   wing_tr=wing_tr, wing_cr=wing_c_root, l_coc=fuselage_co_l, l_cab=fuselage_ca_l,
                                   l_tail=fuselage_ta_l, fus_d=fuselage_diameter, pax=aircraft_n_pax,
                                   geometry_duct=geometry_duct, geometry_pylon=geometry_pylon,
                                   geometry_control=geometry_control, geometry_support=geometry_support,
                                   geometry_nacelle=geometry_nacelle, geometry_propeller=geometry_propeller,
                                   geometry_ht=[], geometry_vt=[], x_duct=x_vector[length], x_wing=x_wing, comp_pe=comp_pe,
                                   propulsor_type=propulsion_type, rpm=RPM, z_duct=z_PE, cv_mode=cv_mode)
        x_cg_vector.append(duuc3.x_cog()[0])

    # -----                     CREATE ATR INSTANCE                    ------ #
    atr: Aircraft = Aircraft(aircraft_type="conventional", conditions=conditions,
                             power=power_condition, bem_input=bem_input, delta_e=delta_e,
                             delta_r=delta_r, reference=reference,
                             wing_span=wing_span, wing_sweep=wing_phi_qc, wing_airfoil=wing_airfoil,
                             wing_tr=wing_tr, wing_cr=wing_c_root, l_coc=fuselage_co_l, l_cab=fuselage_ca_l,
                             l_tail=fuselage_ta_l, fus_d=fuselage_diameter, pax=aircraft_n_pax,
                             geometry_duct=[], geometry_pylon=[], geometry_control=[],
                             geometry_support=[], geometry_nacelle=geometry_nacelle_conv,
                             geometry_propeller=geometry_propeller_conv,
                             geometry_ht=geometry_ht, geometry_vt=geometry_vt, x_duct=0, x_wing=x_wing, comp_pe=[],
                             propulsor_type=propulsion_type, rpm=RPM, z_duct=0, cv_mode=cv_mode)

    # ------                    XCOG EXTRACTION                 ----- #
                 # fus              wing            cg              lemac
    x_cog_duuc = [duuc.x_cog()[2], duuc.x_cog()[3], duuc.x_cog()[0], duuc.x_cog()[1]]
    x_cog_atr = [atr.x_cog()[2], atr.x_cog()[3], atr.x_cog()[0], atr.x_cog()[1]]
    z_cog_duuc = [duuc.x_cog()[4]]
    z_cog_atr = [atr.x_cog()[4]]

    # ------                  VERTICAL TAIL SIZING COUPLED                        ----- #
    eta_h = 0.9
    cd_pe = duuc.empennage.cd_sum()  # drag of one PE during OEI
    cd_wind = 0.06  # drag due to windmilling of the propeller in OEI
    cy_pe = duuc.cy_beta()
    array = np.linspace(0.1, 3.0, 301)

    s_array = []
    s_stab_duuc = []
    fuselage_length = fuselage_ca_l + fuselage_ta_l + fuselage_co_l
    y_center_duuc = np.radians(cant_angle) * (pylon_length + 0.5 * support_length)

    for y in range(len(array)):
        s_array.append(s_control("DUUC", wing_phi_qc, (x_PE - duuc.x_cog()[0]), 2051 * 10 ** 3, eta_h, config.v_approach, 2,
                                 y_center_duuc, cy_duuc=array[y], cd_pe=cd_pe, cd_wind=0.06))
        s_stab_duuc.append(s_stability("DUUC", ref.s_w, duuc.x_cog()[0], fuselage_length, fuselage_diameter,
                           wing_span, (x_PE - duuc.x_cog()[0]), config.v_crit, 0, mach,
                           cl_a_duuc=array[y]))

    surf1 = s_stability("DUUC", ref.s_w, duuc.x_cog()[0], fuselage_length, fuselage_diameter,
                        wing_span, (x_PE - duuc.x_cog()[0]), config.v_crit, 0, mach,
                        cl_a_duuc=duuc.empennage.duct.cl_da())
    surf2 = s_control("DUUC", wing_phi_qc, (x_PE - duuc.x_cog()[0]), 2051 * 10 ** 3, eta_h, config.v_approach, 2,
                      y_center_duuc, cy_duuc=cy_pe, cd_pe=cd_pe, cd_wind=0.06)
    s_vert_req = max(surf1, surf2)

    # ------                  HORIZONTAL TAIL SIZING COUPLED                       ----- #
    cmac = ref.c_mac_w
    a1_atr, b1_atr = slopes("control", "conventional", 0, ref.phi_qc_w, ref.ar_w, fuselage_length,
                            x_cog_atr[3], cmac, eta_h, -0.5, 1.44, 0.597, 0, ref.phi_hc_h, mach,
                            0, ref.phi_hc_h, 0, 0, velocity, ref.s_w, -1.97, 0, 0)
    a2_atr, b2_atr = slopes("stability", "conventional", 3.4, ref.phi_qc_w, ref.ar_w, fuselage_length,
                            x_cog_atr[3], cmac, eta_h, -0.5, 1.44, 0.597, ref.b_w, ref.phi_hc_h, mach, ref.ar_h,
                            ref.phi_hc_h, 0, 0, velocity, ref.s_w, 0, 0, 0)
    z_duuc = z_PE + (pylon_length + 0.5 * support_length) * np.sin(np.radians(cant_angle))

    thrust_duuc = BEM1
    cl_duuc = duuc.empennage.cl_sum()
    cl_a_duuc = duuc.empennage.cl_a()
    cm_h = duuc.empennage.cm_emp()
    a1_duuc, b1_duuc = slopes("control", "DUUC", 0, ref.phi_qc_w, ref.ar_w, fuselage_length,
                              x_cog_duuc[3], cmac, eta_h, 0, 1.44, 0.597, 0, ref.phi_hc_h, mach,
                              0, 0, cm_h, thrust_duuc, velocity, ref.s_w, z_duuc, cl_a_duuc, cl_duuc)
    a2_duuc, b2_duuc = slopes("stability", "DUUC", z_duuc, ref.phi_qc_w, ref.ar_w, fuselage_length,
                              x_cog_duuc[3], cmac, eta_h, 0, 1.44, 0.597, ref.b_w, ref.phi_hc_h, mach,
                              0, 0, 0, thrust_duuc, velocity, ref.s_w, 0, cl_a_duuc, cl_duuc)

    x_cg_duuc = x_cog_duuc[2] - x_cog_duuc[3] - 0.25
    s_w = 61
    static_margin = stm / 100

    def compute_s_hor_req(a2, xcg):
        y_intersection = a2 * xcg + static_margin

        new_sh = y_intersection * s_w

        return new_sh
    s_hor_req = compute_s_hor_req(a2_duuc, x_cg_duuc)

    # -----                 VECTOR EXTRACTION                   ------ #
    w_duuc = [duuc.fuselage.weight(), duuc.wing.weight(), duuc.empennage.duct.weight(),
              duuc.empennage.weight_ps()[0], duuc.empennage.weight_ps()[1],
              duuc.empennage.elevator.weight(), duuc.empennage.propeller.weight_engine(),
              duuc.empennage.nacelle.weight(), duuc.empennage.propeller.weight_fan()]

    w_atr = [atr.fuselage.weight(), atr.wing.weight(), atr.empennage.propeller.weight_engine(),
             atr.empennage.ht_tail.weight(), atr.empennage.vt_tail.weight(), atr.empennage.weight_cv(),
             atr.empennage.weight_nac() / 2, atr.empennage.propeller.weight_engine() / 2,
             atr.empennage.propeller.weight_fan()]

    v_inflow = [velocity, duuc.empennage.duct.inflow_velocity(), duuc.empennage.propeller.inflow_velocity(),
                duuc.empennage.support.inflow_velocity(), duuc.empennage.support.inflow_velocity(),
                duuc.empennage.elevator.inflow_velocity(), duuc.empennage.propeller.u_wake()]

    a_inflow = [alpha, duuc.empennage.duct.inflow_angle()[0], duuc.empennage.propeller.inflow_angle()[0],
                duuc.empennage.support.inflow_angle()[0], duuc.empennage.support.inflow_angle()[0],
                duuc.empennage.elevator.inflow_angle()[0], delta_e / 2]

    pylon_in = [duuc.empennage.pylon.inflow_velocity(), duuc.empennage.pylon.inflow_angle()[0]]

    station = [0, 1, 1.8, 2.2, 2.55, 4.5, 6]
    results = {"Weight": {'w_vector_duuc': w_duuc, 'w_vector_atr': w_atr,
                          "fuselage": [duuc.fuselage.weight(), atr.fuselage.weight()], "wing": [duuc.wing.weight(),
                                                                                                atr.wing.weight()],
                          "empennage": [duuc.empennage.weight(), atr.empennage.weight()]},
               "Inflow": {'v_inflow': v_inflow, 'a_inflow': a_inflow, 'station': station, 'pylon_in': pylon_in},
               "X_cog": {'x_cog_duuc': x_cog_duuc, 'x_cog_atr': x_cog_atr, "z_cog_duuc": z_cog_duuc,
                         "z_cog_atr": z_cog_atr},
               "CD0": {'cd0_duuc': duuc.cd0_empennage(), 'cd0_atr': atr.cd0_empennage()},
               "Vtail": [s_array, s_stab_duuc, 1, 1, 0.95],
               "Htail": {"ATR": [a1_atr, b1_atr, a2_atr], "DUUC": [a1_duuc, b1_duuc, a2_duuc, stm]},
               "Requirements": {"surfaces": [s_vert_req, s_hor_req],
                                "forces_duuc": [duuc.fuselage.lift() + duuc.wing.lift(), duuc.empennage.lift(),
                                                duuc.fuselage.drag() + duuc.wing.drag(), duuc.empennage.drag(),
                                                duuc.thrust()],
                                "forces_atr": [atr.fuselage.lift() + atr.wing.lift(), atr.empennage.lift(),
                                               atr.fuselage.drag() + atr.wing.drag(),
                                               atr.empennage.drag(),
                                               atr.thrust()],
                                "vectors_duuc": [duuc.fuselage.cl() + duuc.wing.cl(), duuc.empennage.cl_sum(),
                                                 duuc.fuselage.cd() + duuc.fuselage.cd(), duuc.empennage.cd_sum(),
                                                 duuc.tc()],
                                "vectors_atr": [atr.fuselage.cl() + atr.wing.cl(), atr.empennage.cl_sum(),
                                                atr.fuselage.cd() + atr.fuselage.cd(), atr.empennage.cd_sum_norm(),
                                                atr.tc()],
                                "w_total": [duuc.weight() * 9.81, atr.weight() * 9.81],
                                "x_cg": [x_cg_vector],
                                "deltas": {"DUUC": [duuc.cm_a(), duuc.cm_de(), duuc.cl_de(), duuc.cl_a(), duuc.cn_dr(),
                                                    duuc.cy_beta(), duuc.cn_beta(), duuc.cy_dr(), duuc.cm_a_PE()],
                                           "ATR": [atr.cm_a(), atr.cm_de(), atr.cl_de(), atr.cl_a(), atr.cn_dr(),
                                                   atr.cy_beta(), atr.cn_beta(), atr.cy_dr(), atr.wing.cm0()]}},
               "Geometry": [duct_profile, pylon_profile, support_profile, cv_airfoil],
               "Pylon": {"Inflow": [duuc.empennage.pylon.inflow_angle()[0], duuc.empennage.pylon.inflow_velocity()],
                         "Cl": [duuc.empennage.pylon.cl()[0], duuc.empennage.pylon.cl()[1]],
                         "Cd": [duuc.empennage.pylon.cd()[0], duuc.empennage.pylon.cd()[1]],
                         "Cm": [duuc.empennage.pylon.cm()[0], duuc.empennage.pylon.cm()[1]]},
               "Support": {"Inflow": [duuc.empennage.support.inflow_angle()[0], duuc.empennage.support.inflow_velocity()],
                           "Cl": [duuc.empennage.support.cl()[0], duuc.empennage.support.cl()[1]],
                           "Cd": [duuc.empennage.support.cd()[0], duuc.empennage.support.cd()[1]],
                           "Cm": [duuc.empennage.support.cm()[0], duuc.empennage.support.cm()[1]]},
               "Duct": {"Inflow": [duuc.empennage.duct.inflow_angle()[0], duuc.empennage.duct.inflow_velocity()],
                        "Cl": [duuc.empennage.duct.cl()[0], duuc.empennage.duct.cl()[1]],
                        "Cd": [duuc.empennage.duct.cd()[0], duuc.empennage.duct.cd()[1]],
                        "Cm": [duuc.empennage.duct.cm()[0], duuc.empennage.duct.cm()[1]]},
               "Nacelle": {"Inflow": [duuc.empennage.nacelle.inflow_angle()[0], duuc.empennage.nacelle.inflow_velocity()],
                           "Cl": [duuc.empennage.nacelle.cl()[0], duuc.empennage.nacelle.cl()[1]],
                           "Cd": [duuc.empennage.nacelle.cd()[0], duuc.empennage.nacelle.cd()[1]],
                           "Cm": [duuc.empennage.nacelle.cm()[0], duuc.empennage.nacelle.cm()[1]],
                           "Vector": [duuc.empennage.nacelle.cl()[2], duuc.empennage.nacelle.cd()[2],
                                      duuc.empennage.nacelle.cm()[2]]},
               "Control": {"Inflow": [duuc.empennage.elevator.inflow_angle()[0], duuc.empennage.elevator.inflow_velocity(),
                                      duuc.empennage.rudder.inflow_angle()[0], duuc.empennage.rudder.inflow_velocity()],
                           "Cl": [duuc.empennage.elevator.cl()[0], duuc.empennage.elevator.cl()[1]],
                           "Cd": [duuc.empennage.elevator.cd()[0], duuc.empennage.elevator.cd()[1],
                                  duuc.empennage.rudder.cd()[0], duuc.empennage.rudder.cd()[1]],
                           "Cm": [duuc.empennage.elevator.cm()[0], duuc.empennage.elevator.cm()[1],
                                  duuc.empennage.rudder.cm()[0], duuc.empennage.rudder.cm()[1]],
                           "Cy": [duuc.empennage.rudder.cl()[0], duuc.empennage.rudder.cl()[1]]},
               "Empennage": {"Vectors": [pe_cl_vector, pe_cd_vector, pe_cm_vector, duuc.empennage.cd_norm_vector(),
                                         duuc.empennage.cl_norm_vector(), duct_vect, pylon_vect, support_vect,
                                         control_vect, nacelle_vect],
                             "Cl": [duuc.empennage.cl_sum()], "Cd": [duuc.empennage.cd_sum()],
                             "Cm": [duuc.empennage.cm_emp()], "Inflow": [duuc.empennage.inflow_angle()]}
               }
    return results


