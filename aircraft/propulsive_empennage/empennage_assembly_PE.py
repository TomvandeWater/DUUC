from analysis_modules.factors import *
from analysis_modules.ISA import air_density_isa
from aircraft.propulsive_empennage.component_properties.pylon import Pylon
from aircraft.propulsive_empennage.component_properties.nacelle import Nacelle
from aircraft.propulsive_empennage.component_properties.support import SupportStrut
from aircraft.propulsive_empennage.component_properties.duct import Duct
from aircraft.propulsive_empennage.component_properties.control_vane import ControlVane
from aircraft.propulsive_empennage.component_properties.propeller import Propeller
import matplotlib.pyplot as plt
import config
import data.atr_reference as ref
from data.read_data import *


class PropulsiveEmpennage:
    def __init__(propemp, rpm: float, alpha: float, power_condition: str, va_inlet: float,
                 n_blades: float, prop_diameter: float, hub_diameter: float,
                 prop_airfoil: str, prop_sweep: float, prop_pitch: float, c_root: float,
                 c_tip: float, duct_diameter: float, duct_chord: float, duct_profile: str,
                 cant_angle: float, pylon_length: float, pylon_chord: float, pylon_profile: str,
                 d_exit: float, nacelle_length: float, nacelle_diameter: float,
                 support_length: float, support_chord: float, support_profile: str, hcv_span: float,
                 hcv_chord: float, control_profile: str, vcv_span: float, vcv_chord: float,
                 propulsor_type: str, v_inf: float, mach: float, ref_area: float, ref_chord: float, ar_wing: float,
                 cl_wing: float, cla_wing: float, bem_input, altitude: float, delta_e: float, delta_r: float):
        super().__init__()
        # operating conditions
        propemp.rpm = rpm
        propemp.ar_wing = ar_wing
        propemp.cl_wing = cl_wing
        propemp.cla_wing = cla_wing
        propemp.alpha = alpha
        propemp.bem_input = bem_input
        propemp.v_inf = v_inf
        propemp.altitude = altitude
        propemp.density = air_density_isa(altitude)[0]
        propemp.delta_e = delta_e
        propemp.delta_r = delta_r

        """ ------------------- determining the inflow angle based on wing downwash -------------------------------- """
        e0 = (2 * propemp.cl_wing) / (np.pi * propemp.ar_wing)
        de_da = 2 * propemp.cla_wing / (np.pi * propemp.ar_wing)
        eta = e0 + de_da * propemp.alpha

        propemp.a_inflow = propemp.alpha - ref.alpha_install_wing - eta + config.ai_duct
        propemp.a_inflow = propemp.alpha

        propemp.power_condition = power_condition
        propemp.va_inlet = area_ratio("0012", duct_chord, duct_diameter/2, 0)[0] * propemp.v_inf
        propemp.cant = cant_angle
        propemp.d_exit = d_exit
        propemp.prop_type = propulsor_type
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
        propemp.ref_chord = ref_chord

        # Initiate propeller class
        propemp.propeller = Propeller(propemp.n_bl, propemp.d_prop, propemp.d_hub,
                                      propemp.prop_airfoil, propemp.s_prop, propemp.p_prop,
                                      propemp.rpm, propemp.power_condition, propemp.va_inlet,
                                      propemp.a_inflow, propemp.ref_area, propemp.c_root,
                                      propemp.c_tip, propemp.v_inf, propemp.prop_type, altitude=propemp.altitude)

        # calculate properties based on input from propeller
        tc_prop = propemp.bem_input[3]
        cn_prop = propemp.propeller.cn()
        va = 2  # """-> should be linked from BEM input """

        # Initiate duct class
        propemp.duct = Duct(propemp.d_duct, propemp.c_duct, propemp.airfoil_duct, propemp.a_inflow,
                            propemp.power_condition, tc_prop, propemp.v_inf,
                            propemp.mach, propemp.ref_area, propemp.ref_chord, propemp.bem_input, va,
                            altitude=propemp.altitude)

        v_ax = 2 * abs(propemp.bem_input[5]) + propemp.v_inf
        v_tan = abs(propemp.bem_input[6])

        angle = np.arctan(v_ax / v_tan)
        v_effective = np.sqrt(v_ax ** 2 + v_tan ** 2)

        v_after_prop = v_effective
        a_after_prop = angle

        # Initiate nacelle class
        propemp.nacelle = Nacelle(propemp.l_nacelle, propemp.d_nacelle, propemp.prop_type,
                                  propemp.power_condition, v_after_prop, propemp.a_inflow, propemp.ref_area,
                                  propemp.v_inf, propemp.mach, a_after_prop, altitude=propemp.altitude)

        # Initiate control vane class for elevator (1 piece)
        propemp.elevator = ControlVane(propemp.b_hcv, propemp.c_hcv, propemp.hcv_profile,
                                       propemp.power_condition, propemp.va_inlet,
                                       propemp.d_exit, propemp.ref_area, propemp.delta_e,
                                       propemp.v_inf, propemp.a_inflow, propemp.mach, v_after_prop, a_after_prop,
                                       altitude=propemp.altitude)

        # Initiate control vane class for rudder (1 piece)
        propemp.rudder = ControlVane(propemp.b_vcv, propemp.c_vcv, propemp.vcv_profile,
                                     propemp.power_condition, propemp.va_inlet, propemp.d_exit,
                                     propemp.ref_area, propemp.delta_r, propemp.v_inf,
                                     propemp.a_inflow, propemp.mach, v_after_prop, a_after_prop,
                                     altitude=propemp.altitude)

        m_supported = ((propemp.rudder.weight() * 4) + propemp.propeller.weight_engine()
                       + propemp.propeller.weight_fan() + propemp.nacelle.weight() + propemp.duct.weight())

        # Initiate support class
        propemp.support = SupportStrut(propemp.l_support, propemp.c_support,
                                       propemp.support_profile, propemp.cant,
                                       propemp.power_condition, v_after_prop, alpha, tc_prop,
                                       cn_prop, propemp.ref_area, propemp.v_inf, propemp.va_inlet, propemp.d_prop,
                                       a_after_prop, m_supported, propemp.ref_chord, altitude=propemp.altitude)

        m_supported2 = ((propemp.rudder.weight() * 4) + propemp.propeller.weight_engine()
                        + propemp.propeller.weight_fan() + propemp.nacelle.weight() + propemp.duct.weight()
                        + propemp.support.weight())

        # Initiate pylon class
        propemp.pylon = Pylon(propemp.l_pylon, propemp.c_pylon, propemp.pylon_profile,
                              propemp.power_condition, propemp.cant, propemp.a_inflow, propemp.ref_area,
                              propemp.v_inf, m_supported2, propemp.ref_chord, altitude=propemp.altitude)
    """ -------------------------------- geometric properties ------------------------------------------------------ """
    def area_wet(self):
        s_wet_duct = self.duct.wetted_area()
        s_wet_pylon = self.pylon.area_wetted()
        s_wet_support = self.support.area_wetted()
        s_wet_control = self.rudder.wet_area()

        s_wet_total = 2 * (s_wet_duct + s_wet_pylon + s_wet_support) + 8 * s_wet_control
        return s_wet_total

    """ -------------------------------- coefficients -------------------------------------------------------------- """
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

    def cd_vector(self):
        cd_duct = self.duct.cdi()[1] + self.duct.cd0()
        cd_pylon = self.pylon.cd()
        cd_nacelle = self.nacelle.cd()
        cd_support = self.support.cd()
        return [cd_duct, cd_pylon, cd_nacelle, cd_support]

    def cd_sum(self):
        cd_sum_pe = sum(self.cd_vector()) * (self.ref_area / self.duct.proj_area()) + sum(self.cd_interference_vector())
        return cd_sum_pe

    def cl_sum(self):
        cl_sum_pe = self.duct.cl()[1] + self.pylon.cl() + self.support.cl()
        # print(f"cl duct: {self.duct.cl()[1]}, cl pylon: {self.pylon.cl()[1]}, cl support: {self.support.cl()[1]}")
        return cl_sum_pe

    def cl_prime(self):
        cl_prime_pe = (self.pylon.cl_prime() + self.duct.cl_prime() + self.support.cl_prime()
                       + self.nacelle.cl_prime() + 2 * self.elevator.cl_prime()
                       + self.propeller.cl_prime())
        return cl_prime_pe

    def cd_prime(self):
        cd_prime_pe = (self.pylon.cd_prime() + self.duct.cd_prime() + self.support.cd_prime()
                       + self.nacelle.cd_prime() + 2 * self.elevator.cd_prime()
                       + 2 * self.rudder.cd_prime() + self.propeller.cd_prime())

        cd_interference = sum(self.cd_interference_vector())

        cd_total = cd_prime_pe + cd_interference
        # print(f"cd prime: {cd_prime_pe}, cd interference: {cd_interference}")
        return cd_total

    def cm_emp(self):
        """ assume the moment to be at the center line of the duct at the c/4 location """
        x_pylon = config.x_pylon
        z_pylon = ((config.duct_diameter / 2) * np.sin(np.radians(config.cant_angle))
                   + (config.pylon_length - config.duct_diameter / 2) / 2 * np.sin(np.radians(config.cant_angle)))
        x_support = config.x_support

        """ properties are all normalized with wing area, inflow velocity and cmac of the wing """
        cm_duct = self.duct.cm()[1]
        cm_support = self.support.cm()[1]
        cm_nac = self.nacelle.cm()[1]
        cm_pylon = self.pylon.cm()[1]

        support = x_support * self.support.cn()[1]
        pylon = x_pylon * self.pylon.cn()[1] + z_pylon * self.pylon.ct()[1]

        mf_sum = -support - pylon
        cm_sum = cm_duct + cm_nac + cm_pylon + cm_support

        propeller_moment = 0  # assume the propeller is located at c/4 and there is no moment in the x-z plane

        if self.power_condition == "off":
            cm_pe = cm_sum + mf_sum
            return cm_pe
        if self.power_condition == "on":
            cm_pe = cm_sum + mf_sum + propeller_moment
            return cm_pe
        else:
            print("Empennage assembly PE.py -> wrong power condition input")
            return -9000000000

    def cl_norm_vector(self):
        cl_norm_duct = self.duct.coefficient_norm()[0]
        cl_norm_pylon = self.pylon.coefficient_norm()[0]
        cl_norm_nacelle = self.nacelle.coefficients_norm()[0]
        cl_norm_support = self.support.coefficient_norm()[0]

        cl_norm_tot = cl_norm_duct + cl_norm_pylon + cl_norm_nacelle + cl_norm_support
        return cl_norm_tot

    def cd_norm_vector(self):
        cd_norm_duct = self.duct.coefficient_norm()[1]
        cd_norm_pylon = self.pylon.coefficient_norm()[1]
        cd_norm_nacelle = self.nacelle.coefficients_norm()[1]
        cd_norm_support = self.support.coefficient_norm()[1]

        cd_norm_tot = cd_norm_duct + cd_norm_pylon + cd_norm_nacelle + cd_norm_support
        return cd_norm_tot

    """ -------------------------------------- Forces -------------------------------------------------------------- """
    def drag(self):
        drag_pe = self.cd_prime() * self.v_inf ** 2 * self.ref_area * 0.5 * self.density
        return drag_pe

    def lift(self):
        lift_pe = self.cd_prime() * self.v_inf ** 2 * self.ref_area * 0.5 * self.density
        return lift_pe

    def thrust(self):
        thrust_pe = 2 * self.propeller.thrust()
        return thrust_pe

    """ -------------------------------------- Weight ------------------------------------------------------------- """

    def weight_ps(self):
        # Safety factors
        n_ult = 1.5  # Ultimate load factor
        k_stoot = 1.5  # Impact factor
        k_misc = 1.25  # Miscellaneous factor

        # Convert cant angle to radians
        c_rad = np.radians(self.cant)

        # Total length of the beam
        l_tot = self.l_support + self.l_pylon  # [m]

        """ Material properties -> AL 6061-T6 """
        sigma_allow = 241 * 1e6  # Allowable stress [Pa] (converted from MPa to N/m^2)
        rho = 2700  # Mass per unit length of the beam [density]

        # Aerodynamic forces on pylon and support
        f_pylon = (0.5 * self.density * self.pylon.inflow_velocity() ** 2 * self.pylon.cl() * self.c_pylon
                   * self.l_pylon)  # [N]
        f_support = (0.5 * self.density * self.support.inflow_velocity() ** 2 * self.support.cl()
                     * self.c_support * self.l_support)  # [N]

        # Aerodynamic moment about the root
        m_aero = (0.5 * f_pylon * self.l_pylon + (0.5 * self.l_support + self.l_pylon) * f_support) * np.cos(
            c_rad)  # [N·m]

        # Weight moment about the root
        m_weight = (self.duct.weight() + self.nacelle.weight() + self.propeller.weight_engine()
                    + self.propeller.weight_fan()) * (self.l_pylon + self.l_support * 0.5) * 9.81 * np.cos(
            c_rad)  # [N·m]

        # Sizing moment (maximum of aerodynamic and weight moments)
        m_sizing = max(m_aero, m_weight)  # [N·m]

        # Determine maximum thickness of the airfoil section
        num_list = [int(digit) for digit in self.pylon_profile]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness percentage
        thickness_pylon = (thickness / 100) * self.c_pylon  # [m]

        num_list2 = [int(digit) for digit in self.support_profile]
        thickness_sup = num_list2[2] * 10 + num_list2[3]  # NACA thickness percentage
        thickness_support = (thickness_sup / 100) * self.c_support  # [m]

        height = max(thickness_pylon, thickness_support) * 0.90  # Reduce maximum thickness by a safety factor [m]

        # Calculate width using bending stress formula
        width = (((6 * m_sizing * n_ult * k_stoot) / sigma_allow) / (height ** 2)) * k_misc  # [m]

        # Total mass of the beam
        mass_tot = height * width * l_tot * rho  # [kg]

        # Split mass into pylon and support sections
        m_pylon = mass_of_section(mass_tot, l_tot, 0, self.l_pylon)  # Mass of pylon section [kg]
        m_support = mass_of_section(mass_tot, l_tot, self.l_pylon, l_tot)  # Mass of support section [kg]

        return m_pylon, m_support

    def weight(self):
        """ determines mass of 1 PE -> output in [kg]"""
        w_pe = (self.duct.weight() + self.propeller.weight_fan() + self.propeller.weight_engine()
                + self.nacelle.weight() + 2 * self.rudder.weight() + 2 * self.elevator.weight() + self.weight_ps()[0]
                + self.weight_ps()[1])
        return w_pe

    """ -------------------------------------- NASA prediction model ----------------------------------------------- """
    def cn_dp(self):
        f1 = 3.10
        f2 = .53
        alpha = np.radians(self.alpha)

        cn_dp_nasa = f1 * np.sin(alpha) * (np.cos(alpha) + f2 * self.yv())
        return cn_dp_nasa

    def ct_dp(self):
        f3 = 1.90
        f4 = .92
        alpha = np.radians(self.alpha)

        ct_dp_nasa = f3 * np.sin(alpha) ** 2 + f4 * (self.yv()) ** 2
        return ct_dp_nasa

    def yv(self):
        alpha = np.radians(self.alpha)
        f4 = .92
        f3 = 1.90
        C_T_dp_guess = 0.168

        A = (config.d_exit / 2) ** 2 * np.pi
        A_p = (config.d_exit / 2) ** 2 * np.pi - (config.nacelle_diameter / 2) ** 2 * np.pi

        a_r = A / A_p

        a_c = -1 * np.cos(alpha)/(a_r * f4 +1)
        b = np.cos(alpha)/(a_r * f4 + 1)
        d = a_r * C_T_dp_guess - (a_r * f3 - 1) * (np.sin(alpha) ** 2)
        c = a_r * f4 + 1

        yv_cal = a_c + np.sqrt(b**2 + d / c)
        return yv_cal

    def cm_dp(self):
        alpha = np.radians(self.alpha)
        f5 = .22
        f6 = 1.5
        f7 = 1.87

        yv = self.yv()

        cm_pe = 4 * f5 * np.sin(alpha) * np.cos(alpha) + (f5 * f6 + f7) * yv * np.sin(alpha)

        return cm_pe

""" Test section"""

if __name__ == "__main__":

    a = np.linspace(0, 20, 41)
    cl = []
    cd = []
    cm = []
    cd2 = []
    cd_intereference = []
    ct_dp = []
    cn_dp = []
    cm_dp = []

    exp2_ref_area = 3.51  # 2.57 * 1.2246

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
    cl_exp2_n = []
    cd_exp2 = [-0.042,-0.028,-0.014,0,0.014,0.028,0.042,0.056]
    cd_exp2_n = []

    cd_exp_v11 = [0.40,0.33,0.30,0.296,0.304,0.312,0.32,0.328,0.336,0.344,0.352,0.36,0.368,0.376,0.384,0.392,0.4,0.408,
                  0.416,0.424,0.432,0.44,0.448,0.456,0.464,0.472,0.48,0.488,0.496,0.504,0.512,0.52,0.528,0.536,0.544,
                  0.552,0.56,0.568,0.576,0.584,0.592,0.6,0.608,0.616,0.624,0.632,0.64,0.648,0.656,0.664,0.672,0.68,
                  0.688,0.696,0.704,0.712,0.72,0.728,0.736,0.744]
    cl_exp_v11 = [-0.34663,-0.14754,0.05553,0.11747,0.17941,0.24134,0.30328,0.33325,0.36322,0.39319,0.42316,0.46023,
                  0.49812,0.536,0.57389,0.61178,0.64967,0.68755,0.7117,0.73549,0.75928,0.78307,0.80686,0.83065,0.85444,
                  0.87823,0.90202,0.92581,0.9496,0.96281,0.97545,0.98809,1.00072,1.01336,1.026,1.03864,1.05128,1.06392,
                  1.07655,1.08919,1.09075,1.09139,1.09203,1.09267,1.09331,1.09395,1.09459,1.09523,1.09587,1.09651,
                  1.09715,1.09779,1.09719,1.08577,1.07436,1.06294,1.05152,1.04011,1.02869,1.01727]


    cl_7foot2 = [0.35,1.436,1.844,2.18,2.435,2.588,2.723,2.842,2.941]
    cd_7foot2= [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6]

    cl_7foot =[0.03485,0.53909,0.75034,0.90017,1.0234,1.06834,1.11329,1.15823,1.14826,0.99863]
    cd_7foot = [-5.6,-4.2,-2.8,-1.4,0,1.4,2.8,4.2,5.6,7]

    # duct reference values for gui
    a_ref_duct = np.linspace(0, 15, 16)
    cl_ref_duct = [0, 0.1, 0.203, 0.291, 0.357, 0.433, 0.511, 0.577, 0.702, 0.768, 0.844, 0.935, 1.067, 1.15, 1.2, 1.3]
    a_exp_duct = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
             29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
    cl_exp_duct = [0.116, 0.198, 0.284, 0.369, 0.454, 0.54, 0.629, 0.719, 0.808, 0.898, 0.988, 1.083, 1.181, 1.28, 1.379,
              1.478, 1.577, 1.677, 1.776, 1.876, 1.976, 2.073, 2.17, 2.266, 2.363, 2.459, 2.446, 2.363, 2.281, 2.198,
              2.115, 2.093, 2.073, 2.053, 2.033, 2.012, 1.965, 1.917, 1.868, 1.82, 1.771, 1.725, 1.678, 1.632, 1.586,
              1.539, 1.492, 1.445, 1.398, 1.351, 1.304]
    cd_ref_duct = [0.027, 0.028, 0.03, 0.033, 0.038, 0.041, 0.047, 0.054, 0.061, 0.071, 0.08, 0.087, 0.096, 0.109,
              0.125, 0.129]

    results_duct = []
    results_pylon = []
    results_support = []
    results_nacelle = []
    results_control_vanes = []
    reference_results_pylon = []
    reference_results_support = []
    reference_results_duct = []
    reference_results_duct2 = []

    cl_7foot_n = []
    cd_7foot_n = []

    cl_7foot2_n = []
    cd_7foot2_n = []

    cd_vikesh_n = []
    cl_vikesh_n = []

    cl_vikesh = [0.44262, -0.39344, 1.09016, 1.0082]
    cd_vikesh = [0.35064, 0.3928, 0.68689, 0.7455]

    for i in range(len(a)):
        PE = PropulsiveEmpennage(rpm=config.rpm,
                                 alpha=a[i],
                                 power_condition=config.power_condition,
                                 va_inlet=723,
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
                                 ref_area=ref.s_w,
                                 ref_chord=2.2345,
                                 ar_wing=12,
                                 cl_wing=1.44,
                                 cla_wing=5.89,
                                 bem_input=[41420.85603250924, 26482.06279555917, -1.4475057750305072e-13, 0.8892292886261024, 0.3297344029147765, 0.3201225439053968, 5, 10],
                                 delta_e=0, delta_r=0, altitude=7000)
        cd_intereference.append(PE.cd_interference_vector())
        cl.append(PE.cl_norm_vector())
        cd.append(PE.cd_norm_vector())
        cd2.append(PE.cd_sum())
        cm.append(PE.cm_emp())
        cm_dp.append(-PE.cm_dp())
        cn_dp.append(PE.cn_dp())

        polar_pylon = airfoil_polar(f"pylon0012.txt", float(a[i]))
        cd_val_pylon = float(polar_pylon[1])
        cl_val_pylon = float(polar_pylon[0])
        cm_val_pylon = float(polar_pylon[2])

        polar_support = airfoil_polar(f"support0012.txt", float(a[i]))
        cd_val_support = float(polar_support[1])
        cl_val_support = float(polar_support[0])
        cm_val_support = float(polar_support[2])

        results_duct.append([a[i], PE.duct.cl()[0], PE.duct.cd(), PE.duct.cm()[0]])
        results_pylon.append([a[i], PE.pylon.cl(), PE.pylon.cd(), PE.pylon.cm()[0]])
        results_support.append([a[i], PE.support.cl(), PE.support.cd(), PE.support.cm()[0]])
        results_nacelle.append([a[i], PE.nacelle.cl(), PE.nacelle.cd(), PE.nacelle.cm()[0]])
        results_control_vanes.append([a[i], PE.elevator.cl(), PE.elevator.cd(), PE.elevator.cm()])

        reference_results_pylon.append([a[i], cl_val_pylon, cd_val_pylon, cm_val_pylon])
        reference_results_support.append([a[i], cl_val_support, cd_val_support, cm_val_support])
        # for l in range(len(cd_exp_v11)):
        #    cl_exp_v2.append(2 * cl_exp_v11[l] * exp2_ref_area / PE.ref_area)
        #    cd_exp_v2.append(2 * cd_exp_v11[l] * exp2_ref_area / PE.ref_area)

        for u in range(len(a_ref_duct)):
            reference_results_duct.append([a_ref_duct[u], cl_ref_duct[u], cd_ref_duct[u], None])

        for o in range(len(a_exp_duct)):
            reference_results_duct2.append([a_exp_duct[o], cl_exp_duct[o], None, None])

        for k in range(len(cl_7foot)):
            cl_7foot_n.append(cl_7foot[k] * (2.985 / 61))
            cd_7foot_n.append(cd_7foot[k] * (2.985 / 61))

        for n in range(len(cl_7foot2)):
            cl_7foot2_n.append(cl_7foot2[n] * (2.985 / 61))
            cd_7foot2_n.append(cd_7foot2[n] * (2.985 / 61))

        for j in range(len(cl_vikesh)):
            cl_vikesh_n.append(cl_vikesh[j] / 3.75)
            cd_vikesh_n.append(cd_vikesh[j] / 3.75)

    alpha_avl, cl_avl, clff_avl, cd_avl, cdin_avl, cdff_avl = read_avl_output("AVL_RW_aeroproperties.txt")

    results = np.array(cd_intereference)

    PE.weight_ps()
"""
    create_gui(results_duct, results_pylon, results_nacelle, results_support, results_control_vanes,
               reference_results_duct, reference_results_pylon, None,
               reference_results_support, None, reference_results_duct2,
               None, None, None, None) """

"""
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
    #plt.plot(a_exp1, cl_exp1, label=r'Experimental 2', color="tab:red", marker='o')
    #plt.plot(a_exp_v, cl_exp_v, label=r'Experimental 3 - power off', color="tab:purple", marker='o')
    plt.plot(alpha_avl, cl_avl, label=r'AVL', color='tab:orange', marker='x')
    # plt.plot(alpha_avl, clff_avl, label=r'AVL', color='tab:orange', marker='x')
    # plt.plot(a, cn_dp, label="nasa", color="orange")
    plt.xlim([0, 20])
    plt.ylim([0, 0.5])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $\alpha$ - Propulsive Empennage')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Model', color="tab:blue")
    #plt.plot(alpha_avl, cd_avl, label=r'AVL', color='tab:orange', marker='x')
    #plt.plot(alpha_avl, cdin_avl, label=r'AVL - ind', color='tab:orange', marker='x')
    plt.plot(alpha_avl, cdff_avl, label=r'AVL', color='tab:orange', marker='x')
    # plt.plot(a, ct_dp, label="nasa", color="orange")
    # plt.plot(a_ref, cd_ref, label=r'Experimental', color="tab:red", marker='o')
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $\alpha$ - Propulsive Empennage')
    plt.legend()
    plt.grid(True)

    plt.figure('CM - alpha')
    plt.plot(a, cm, label=r'Model', color="tab:blue")
    # plt.plot(a, cm_dp, label="nasa", color="orange")
    # plt.plot(a_ref, cd_ref, label=r'Experimental', color="tab:red", marker='o')
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{M}$ [-]')
    plt.title(r'Propulsive Empennage')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CM')
    plt.plot(cm, cl, label=r'Model', color="tab:blue")
    # plt.plot(cm_ref, cl_cm_ref, label="nasa", color="orange")
    plt.xlabel(r'$C_{M}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Propulsive Empennage')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD')
    plt.plot(cd, cl, label=r'Model', color="tab:blue")
    # plt.scatter(cd_7foot2_n, cl_7foot2_n, label=r'Experimental 2', color="tab:green", marker='o')
    plt.plot(cd_avl, cl_avl, label=r'AVL', color='tab:orange', marker='x')
    # plt.plot(ct_dp, cn_dp, label="nasa", color="orange")
    # plt.scatter(cd_7foot_n, cl_7foot_n, label=r'Experimental 2', color="tab:red", marker='o', linestyle="dashed")
    # plt.scatter(cd_vikesh_n, cl_vikesh_n, label=r'Experimental 3', color="tab:purple", marker='o')
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    # plt.xlim([0, 0.04])
    # plt.ylim([-0.1, 0.3])
    plt.title(r'$C_{L}$ vs. $C_{D}$ - Propulsive Empennage')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.show() """


