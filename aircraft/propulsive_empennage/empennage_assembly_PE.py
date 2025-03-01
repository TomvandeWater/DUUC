import flow_conditions
from aircraft.propulsive_empennage.component_properties.pylon import Pylon
from aircraft.propulsive_empennage.component_properties.nacelle import Nacelle
from aircraft.propulsive_empennage.component_properties.support import SupportStrut
from aircraft.propulsive_empennage.component_properties.duct import Duct
from aircraft.propulsive_empennage.component_properties.control_vane import ControlVane
from aircraft.propulsive_empennage.component_properties.propeller import Propeller
import numpy as np
import matplotlib.pyplot as plt
import config
import data.atr_reference as ref


class PropulsiveEmpennage:
    def __init__(propemp, rpm: float, alpha: float, power_condition: str, va_inlet: float,
                 re_duct: float, n_blades: float, prop_diameter: float, hub_diameter: float,
                 prop_airfoil: str, prop_sweep: float, prop_pitch: float, c_root: float,
                 c_tip: float, duct_diameter: float, duct_chord: float, duct_profile: str,
                 cant_angle: float, pylon_length: float, pylon_chord: float, pylon_profile: str,
                 d_exit: float, nacelle_length: float, nacelle_diameter: float,
                 support_length: float, support_chord: float, support_profile: str, hcv_span: float,
                 hcv_chord: float, control_profile: str, vcv_span: float, vcv_chord: float,
                 propulsor_type: str, v_inf: float, mach: float, ref_area: float):
        super().__init__()
        # operating conditions
        propemp.rpm = rpm
        propemp.alpha = alpha
        propemp.power_condition = power_condition
        propemp.va_inlet = va_inlet
        propemp.re = re_duct
        propemp.cant = cant_angle
        propemp.d_exit = d_exit
        propemp.prop_type = propulsor_type
        propemp.v_inf = v_inf
        propemp.mach = mach

        # propeller properties
        propemp.n_bl = n_blades
        propemp.d_prop = prop_diameter
        propemp.d_hub = hub_diameter
        propemp.prop_airfoil = prop_airfoil
        propemp.s_prop = prop_sweep
        propemp.p_prop = prop_pitch
        propemp.c_root = c_root
        propemp.c_tip = c_tip

        # duct properties
        propemp.d_duct = duct_diameter
        propemp.c_duct = duct_chord
        propemp.airfoil_duct = duct_profile

        # pylon properties
        propemp.l_pylon = pylon_length
        propemp.c_pylon = pylon_chord
        propemp.pylon_profile = pylon_profile

        # nacelle properties
        propemp.l_nacelle = nacelle_length
        propemp.d_nacelle = nacelle_diameter

        # support properties
        propemp.l_support = support_length
        propemp.c_support = support_chord
        propemp.support_profile = support_profile

        # horizontal control vane properties
        propemp.b_hcv = hcv_span
        propemp.c_hcv = hcv_chord
        propemp.hcv_profile = control_profile

        # vertical control vane properties
        propemp.b_vcv = vcv_span
        propemp.c_vcv = vcv_chord
        propemp.vcv_profile = control_profile

        propemp.ref_area = ref_area

        # Initiate propeller class
        propemp.propeller = Propeller(propemp.n_bl, propemp.d_prop, propemp.d_hub,
                                      propemp.prop_airfoil, propemp.s_prop, propemp.p_prop,
                                      propemp.rpm, propemp.power_condition, propemp.va_inlet,
                                      propemp.alpha, propemp.re, propemp.ref_area, propemp.c_root,
                                      propemp.c_tip, propemp.v_inf)

        # calculate properties based on input from propeller
        u1 = propemp.propeller.u1()
        # print(f"u1 = {u1}")
        tc_prop = propemp.propeller.tc()
        # print(f"tc_prop = {tc_prop}")
        v_prop = propemp.propeller.inflow_velocity()
        cn_prop = propemp.propeller.cn()
        # print(f"cn_prop = {cn_prop}")
        u_mom = 0.5 * (v_prop + u1)
        # print(f"u_mom: {u_mom}")

        # Initiate duct class
        propemp.duct = Duct(propemp.d_duct, propemp.c_duct, propemp.airfoil_duct, propemp.alpha,
                            propemp.re, propemp.power_condition, u_mom, tc_prop, propemp.v_inf,
                            propemp.mach, propemp.ref_area)

        # Initiate nacelle class
        propemp.nacelle = Nacelle(propemp.l_nacelle, propemp.d_nacelle, propemp.prop_type,
                                  propemp.re, propemp.power_condition, u_mom, alpha, propemp.ref_area,
                                  propemp.v_inf, propemp.mach)

        # v_after_prop = np.sqrt(propemp.propeller.vnorm() ** 2 + propemp.propeller.vax() ** 2)
        v_after_prop = 100
        # a_after_prop = np.arctan(propemp.propeller.vnorm()/propemp.propeller.vax())
        a_after_prop = 10

        # Initiate control vane class for elevator (1 piece)
        propemp.elevator = ControlVane(propemp.b_hcv, propemp.c_hcv, propemp.hcv_profile,
                                       propemp.power_condition, propemp.va_inlet,
                                       propemp.d_exit, u1, propemp.ref_area, flow_conditions.delta_e,
                                       propemp.re, propemp.v_inf, propemp.alpha, propemp.mach, u_mom)

        # Initiate control vane class for rudder (1 piece)
        propemp.rudder = ControlVane(propemp.b_vcv, propemp.c_vcv, propemp.vcv_profile,
                                     propemp.power_condition, propemp.va_inlet, propemp.d_exit, u1,
                                     propemp.ref_area, flow_conditions.delta_r, propemp.re, propemp.v_inf,
                                     propemp.alpha, propemp.mach, u_mom)

        m_supported = ((propemp.rudder.weight() * 0.25 / 2) + propemp.propeller.weight_engine()
                       + propemp.propeller.weight_fan() + propemp.nacelle.weight() + propemp.duct.weight())

        # Initiate support class
        propemp.support = SupportStrut(propemp.l_support, propemp.c_support,
                                       propemp.support_profile, propemp.cant,
                                       propemp.power_condition, v_after_prop, u_mom, alpha, tc_prop,
                                       cn_prop, propemp.ref_area, propemp.v_inf, propemp.va_inlet, propemp.d_prop,
                                       a_after_prop, m_supported)

        # Initiate pylon class
        propemp.pylon = Pylon(propemp.l_pylon, propemp.c_pylon, propemp.pylon_profile,
                              propemp.power_condition, propemp.cant, propemp.alpha, propemp.ref_area,
                              propemp.v_inf, m_supported)

    """ -------------------------------- coefficients ------------------------------------------------------------- """
    def cd0_vector(self):
        cd0_duct = self.duct.cd0()
        cd0_pylon = self.pylon.cd0()
        cd0_nacelle = self.nacelle.cd0()
        cd0_support = self.support.cd0()
        cd0_control = self.elevator.cd0()
        return [cd0_duct, cd0_pylon, cd0_nacelle, cd0_control, cd0_control, cd0_support]

    def cd_interference_vector(self):
        cd_interference_pylon = self.pylon.cd_interference() * 2
        cd_interference_support = self.support.cd_interference() * 2
        cd_interference_control = self.rudder.cd_interference() * 8
        cd_interference_propeller = self.propeller.cd_interference() * 2

        return [cd_interference_pylon, cd_interference_support, cd_interference_control, cd_interference_propeller]

    def cl_prime(self):
        cl_prime_pe = (self.pylon.cl_prime() + self.duct.cl_prime() + self.support.cl_prime()
                       + self.nacelle.cl_prime() + 2 * self.elevator.cl_prime()
                       + self.propeller.cl_prime())
        return cl_prime_pe

    def cd_prime(self):
        cd_prime_pe = (self.pylon.cd_prime() + self.duct.cd_prime() + self.support.cd_prime()
                       + self.nacelle.cd_prime() + 2 * self.elevator.cd_prime()
                       + 2 * self.rudder.cd_prime() + self.propeller.cd_prime())

        cd_interference = (self.pylon.cd_interference() + self.support.cd_interference()
                           + self.propeller.cd_interference() + 2 * self.rudder.cd_interference()
                           + 2 * self.elevator.cd_interference())

        cd_total = cd_prime_pe + cd_interference
        # print(f"cd prime: {cd_prime_pe}, cd interference: {cd_interference}")
        return cd_total

    def cm_emp(self):
        l_duct = 0.25 * self.duct.duct_chord
        l_support = 0.45 * self.duct.duct_chord + (0.25 * self.support.support_chord)
        l_pylon = 0.45 * self.duct.duct_chord + (0.25 * self.pylon.pylon_chord)

        cm_duct = self.duct.cm()
        cm_support = self.support.cm()
        cm_nac = self.nacelle.cm()
        cm_pylon = self.pylon.cm()
        l_duct_dist = self.cl_prime() * l_duct
        l_sup_dist = self.support.cl_prime() * l_support
        l_pyl_dist = self.pylon.cl_prime() * l_pylon

        cm_pe = cm_duct + cm_support + cm_nac + cm_pylon - l_duct_dist - l_sup_dist - l_pyl_dist

        return cm_pe

    """ -------------------------------------- Forces -------------------------------------------------------------- """
    def drag(self):
        drag_pe = self.cd_prime() * self.v_inf ** 2 * self.ref_area * 0.5 * flow_conditions.rho
        return drag_pe

    def lift(self):
        lift_pe = self.cd_prime() * self.v_inf ** 2 * self.ref_area * 0.5 * flow_conditions.rho
        return lift_pe

    def thrust(self):
        thrust_pe = 2 * self.propeller.thrust()
        return thrust_pe

    """ -------------------------------------- Weight ------------------------------------------------------------- """
    def weight(self):
        """ determines mass of 1 PE -> output in [kg]"""
        w_pe = (self.duct.weight() + self.propeller.weight_fan() + self.propeller.weight_engine()
                + self.nacelle.weight() + 2 * self.rudder.weight() + 2 * self.elevator.weight() + self.pylon.weight()
                + self.support.weight())
        return w_pe


""" Test section"""

if __name__ == "__main__":

    a = np.linspace(0, 20, 41)
    cl = []
    cd = []
    cd_intereference = []

    cl_exp1 = [0.095,0.182,0.267,0.347,0.426,0.522,0.62,0.707,0.777,0.846,0.907,0.965,1.024,1.078,1.131,1.158,1.151,
               1.145,1.129,1.108,1.074,1.04,1.006]
    a_exp1 = [4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92]

    a_exp_v = [-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
    cl_exp_11 = [-0.38,-0.299,-0.218,-0.137,-0.056,0.025,0.107,0.188,0.269,0.35,0.431,0.499,0.565,0.631,0.697,0.763,
                 0.828,0.894,0.96,1.026,1.092,1.077,1.06,1.044,1.027]

    cl_exp_v = []
    for i in range(len(cl_exp_11)):
        cl_exp_v.append((cl_exp_11[i]/10) * 2)

    cl_exp2 = [0.333,0.678,0.881,0.962,1.032,1.092,1.152,1.126]
    cd_exp2 = [-0.042,-0.028,-0.014,0,0.014,0.028,0.042,0.056]

    cd_exp_v11 = [0.40,0.33,0.30,0.296,0.304,0.312,0.32,0.328,0.336,0.344,0.352,0.36,0.368,0.376,0.384,0.392,0.4,0.408,
                  0.416,0.424,0.432,0.44,0.448,0.456,0.464,0.472,0.48,0.488,0.496,0.504,0.512,0.52,0.528,0.536,0.544,
                  0.552,0.56,0.568,0.576,0.584,0.592,0.6,0.608,0.616,0.624,0.632,0.64,0.648,0.656,0.664,0.672,0.68,
                  0.688,0.696,0.704,0.712,0.72,0.728,0.736,0.744]
    cl_exp_v11 = [-0.34663,-0.14754,0.05553,0.11747,0.17941,0.24134,0.30328,0.33325,0.36322,0.39319,0.42316,0.46023,
                  0.49812,0.536,0.57389,0.61178,0.64967,0.68755,0.7117,0.73549,0.75928,0.78307,0.80686,0.83065,0.85444,
                  0.87823,0.90202,0.92581,0.9496,0.96281,0.97545,0.98809,1.00072,1.01336,1.026,1.03864,1.05128,1.06392,
                  1.07655,1.08919,1.09075,1.09139,1.09203,1.09267,1.09331,1.09395,1.09459,1.09523,1.09587,1.09651,
                  1.09715,1.09779,1.09719,1.08577,1.07436,1.06294,1.05152,1.04011,1.02869,1.01727]

    cd_exp_v2 = []
    cl_exp_v2 = []

    for i in range(len(cd_exp_v11)):
        cl_exp_v2.append(2 * (cl_exp_v11[i]/10))
        cd_exp_v2.append(2 * (cd_exp_v11[i] / 10))

    for i in range(len(a)):
        PE = PropulsiveEmpennage(rpm=config.rpm,
                                 alpha=a[i],
                                 power_condition=config.power_condition,
                                 va_inlet=723,
                                 re_duct=84222274,
                                 n_blades=config.n_blades,
                                 prop_diameter=config.duct_diameter * 0.95,
                                 hub_diameter=config.hub_diameter,
                                 prop_airfoil=config.prop_airfoil,
                                 prop_sweep=config.propeller_sweep,
                                 prop_pitch=config.propeller_pitch,
                                 c_root=config.c_root,
                                 c_tip=config.c_tip,
                                 duct_diameter=config.duct_diameter,
                                 duct_chord=config.duct_chord,
                                 duct_profile=config.duct_airfoil,
                                 cant_angle=config.cant_angle,
                                 pylon_length=config.pylon_length,
                                 pylon_chord=config.pylon_chord,
                                 pylon_profile=config.pylon_airfoil,
                                 d_exit=config.d_exit,
                                 nacelle_length=config.nacelle_length,
                                 nacelle_diameter=config.nacelle_diameter,
                                 support_length=config.support_length,
                                 support_chord=config.support_chord,
                                 support_profile=config.support_airfoil,
                                 hcv_span=config.control_vane_length,
                                 hcv_chord=config.control_vane_chord,
                                 control_profile=config.control_vanes_airfoil,
                                 vcv_span=config.control_vane_length,
                                 vcv_chord=config.control_vane_chord,
                                 propulsor_type=config.propulsor_type,
                                 v_inf=128,
                                 mach=0.44,
                                 ref_area=ref.s_w)
        cd_intereference.append(PE.cd_interference_vector())
        cl.append(PE.cl_prime())
        cd.append(PE.cd_prime())

    results = np.array(cd_intereference)

    plt.figure('Interference drag')
    plt.plot(a, results[:, 0], label=r'Pylon')
    plt.plot(a, results[:, 1], label=r'Support')
    plt.plot(a, results[:, 2], label=r'Control')
    plt.plot(a, results[:, 3], label=r'Propeller')
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D_{interference}}$ [-]')
    plt.title(r'$C_{D_{interference}}$ vs. $\alpha$ - Propulsive Empennage')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - alpha')
    plt.plot(a, cl, label=r'Model', color="tab:blue")
    plt.plot(a_exp1, cl_exp1, label=r'Experimental 2', color="tab:red", marker='o')
    plt.plot(a_exp_v, cl_exp_v, label=r'Experimental 3', color="tab:purple", marker='o')
    plt.xlim([0, 20])
    plt.ylim([0, 0.5])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $\alpha$ - Propulsive Empennage')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Model', color="tab:blue")
    # plt.plot(a_ref, cd_ref, label=r'Experimental', color="tab:red", marker='o')
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $\alpha$ - Propulsive Empennage')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD')
    plt.plot(cd, cl, label=r'Model', color="tab:blue")
    plt.plot(cd_exp2, cl_exp2, label=r'Experimental 2', color="tab:red", marker='o')
    plt.plot(cd_exp_v2, cl_exp_v2, label=r'Experimental 3', color="tab:purple", marker='o')
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    # plt.xlim([0, 0.25])
    # plt.ylim([-0.1, 1.6])
    plt.title(r'$C_{L}$ vs. $C_{D}$ - Propulsive Empennage')
    plt.legend()
    plt.grid(True)
    plt.show()
