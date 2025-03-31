from aircraft.aircraft_assembly import Aircraft
from analysis_modules.aerodynamic import speed_of_sound
from analysis_modules.vtail_sizing import *


def calculation_manager(parameters):
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

    bem_input = [BEM1, BEM2, BEM3, BEM4, BEM5, BEM6, BEM7, BEM8]
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

    fuselage_diameter = parameters.get('fuselage_diameter')
    fuselage_co_l = parameters.get('fuselage_co_l')
    fuselage_ca_l = parameters.get('fuselage_ca_l')
    fuselage_ta_l = parameters.get('fuselage_ta_l')

    aircraft_n_pax = parameters.get('aircraft_n_pax')
    x_PE = parameters.get('x_PE')
    y_PE = parameters.get('y_PE')
    z_PE = parameters.get('z_PE')

    nacelle_length = parameters.get('nacelle_length')
    nacelle_diameter = parameters.get('nacelle_diameter')

    hcv_span = parameters.get('hcv_span')
    hcv_chord = parameters.get('hcv_chord')

    cv_airfoil = parameters.get('cv_airfoil')

    duuc: Aircraft = Aircraft(aircraft_type="DUUC", alpha=alpha, v_inf=velocity, mach=mach,
                              power=power_condition, bem_input=bem_input, delta_e=delta_e,
                              delta_r=delta_r, altitude=altitude,
                              wing_span=wing_span, wing_sweep=wing_phi_qc, wing_airfoil=wing_airfoil,
                              wing_tr=wing_tr, wing_cr=wing_c_root, l_coc=fuselage_co_l, l_cab=fuselage_ca_l,
                              l_tail=fuselage_ta_l, fus_d=fuselage_diameter, pax=aircraft_n_pax,
                              d_duct=duct_diameter, c_duct=duct_chord, duct_airfoil=duct_profile,
                              b_pylon=pylon_length, c_pylon=pylon_chord, pylon_airfoil=pylon_profile,
                              cant_angle=cant_angle, b_support=support_length, c_support=support_chord,
                              support_airfoil=support_profile, cv_airfoil=cv_airfoil, n_blades=num_blades,
                              l_nacelle=nacelle_length, d_nacelle=nacelle_diameter, l_cv=hcv_span,
                              c_cv=hcv_chord, propulsor_type=propulsion_type, prop_airfoil=propeller_airfoil,
                              d_hub=hub_diameter, rpm=RPM, prop_c_root=propeller_c_root, prop_c_tip=propeller_c_tip,
                              prop_diameter=propeller_diameter)

    atr: Aircraft = Aircraft(aircraft_type="conventional", alpha=alpha, v_inf=velocity, mach=mach,
                             power=power_condition, bem_input=bem_input, delta_e=delta_e,
                             delta_r=delta_r, altitude=altitude,
                             wing_span=ref.b_w, wing_sweep=ref.phi_qc_w, wing_airfoil=ref.wing_airfoil,
                             wing_tr=ref.tr_w, wing_cr=ref.c_root_w, l_coc=ref.l_cockpit, l_cab=ref.l_cab,
                             l_tail=ref.l_tail, fus_d=ref.diameter_fuselage, pax=aircraft_n_pax,
                             d_duct=duct_diameter, c_duct=duct_chord, duct_airfoil=duct_profile,
                             b_pylon=pylon_length, c_pylon=pylon_chord, pylon_airfoil=pylon_profile,
                             cant_angle=cant_angle, b_support=support_length, c_support=support_chord,
                             support_airfoil=support_profile, cv_airfoil=cv_airfoil, n_blades=num_blades,
                             l_nacelle=ref.l_nacelle, d_nacelle=ref.d_nacelle, l_cv=hcv_span,
                             c_cv=hcv_chord, propulsor_type=propulsion_type, prop_airfoil=propeller_airfoil,
                             d_hub=hub_diameter, rpm=RPM, prop_c_root=propeller_c_root, prop_c_tip=propeller_c_tip,
                             prop_diameter=propeller_diameter)

    # ------                  Vertical Tail sizing coupled                        ----- #
    array = np.linspace(0.1, 3.0, 301)

    s_array = []
    s_array_atr = []
    s_stab_duuc = []
    s_stab_atr = []
    fuselage_length = fuselage_ca_l + fuselage_ta_l + fuselage_co_l

    for y in range(len(array)):
        s_array.append(s_control("DUUC", wing_phi_qc, x_PE, 2051 * 10 ** 3, 0.73, 60.4, 2, 2.8,
                                 cy_duuc=array[y], cd_pe=0.00568, cd_wind=0.06))
        s_array_atr.append(s_control("conventional", wing_phi_qc, x_PE, 2051 * 10 ** 3, 0.73, 60.4, 2, 3.608,
                                     cy_atr=array[y]))
        s_stab_duuc.append(s_stability("DUUC", ref.s_w, duuc.x_cog()[0], fuselage_length, fuselage_diameter,
                                       ref.b_w, x_PE, 141, ref.ar_v, 0.31, mach, cy_duuc=array[y]))
        s_stab_atr.append(s_stability("conventional", ref.s_w, atr.x_cog()[0], fuselage_length, fuselage_diameter,
                                      ref.b_w, x_PE, 141, ref.ar_v, 0.31, mach, cy_duuc=array[y]))

    cyd = 0.975
    a = 13.788741676315057
    b = s_control("DUUC", wing_phi_qc, x_PE, 2051 * 10 ** 3, 0.73, 60.4, 2, 2.8,
                  cy_duuc=cyd, cd_pe=0.00568, cd_wind=0.06)
    s_req = max(a, b)


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
                duuc.empennage.elevator.inflow_velocity(), duuc.empennage.propeller.u1(), velocity]

    a_inflow = [alpha, duuc.empennage.duct.inflow_angle(), duuc.empennage.propeller.inflow_angle(),
                duuc.empennage.support.inflow_angle(), duuc.empennage.support.inflow_angle(),
                duuc.empennage.elevator.inflow_angle(), delta_e / 2, alpha]

    station = [0, 1, 1.8, 2.2, 2.55, 4.5, 5.5, 6]
                 # fus              wing            cg              lemac
    x_cog_duuc = [duuc.x_cog()[2], duuc.x_cog()[3], duuc.x_cog()[0], duuc.x_cog()[1]]
    x_cog_atr = [atr.x_cog()[2], atr.x_cog()[3], atr.x_cog()[0], atr.x_cog()[1]]

    results = {"Weight": {'w_vector_duuc': w_duuc, 'w_vector_atr': w_atr},
               "Inflow": {'v_inflow': v_inflow, 'a_inflow': a_inflow, 'station': station},
               "X_cog": {'x_cog_duuc': x_cog_duuc, 'x_cog_atr': x_cog_atr},
               "CD0": {'cd0_duuc': duuc.cd0_empennage(), 'cd0_atr': atr.cd0_empennage()},
               "Vtail": [s_array, s_array_atr, s_stab_duuc, s_stab_atr, a, b, cyd],
               "requirements": [s_req], "Geometry": [0]}
    return results


