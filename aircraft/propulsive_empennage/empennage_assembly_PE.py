from analysis_modules.factors import *
from analysis_modules.ISA import air_density_isa
from aircraft.propulsive_empennage.component_properties.pylon import Pylon
from aircraft.propulsive_empennage.component_properties.nacelle import Nacelle
from aircraft.propulsive_empennage.component_properties.support import SupportStrut
from aircraft.propulsive_empennage.component_properties.duct import Duct
from aircraft.propulsive_empennage.component_properties.control_vane import ControlVane
from aircraft.propulsive_empennage.component_properties.propeller import Propeller
import data.experiment_reference_5annular_airfoil as ref5r
import data.experiment_reference_7foot as ref7f
import data.experiment_reference_linear_regression as reflr
from scipy.interpolate import interp1d
import data.experiment_avl as refavl
import config
import data.atr_reference as ref
from data.read_data import *


class PropulsiveEmpennage:
    def __init__(propemp, conditions, reference, geometry_duct, geometry_pylon, geometry_support, geometry_nacelle,
                 geometry_control, geometry_propeller, rpm: float, ar_wing: float, cl_wing: float, cla_wing: float,
                 bem_input, delta_e: float, delta_r: float, propulsor_type: str, power_condition: str, va_inlet: float,
                 d_exit: float, comp_pe):
        super().__init__()
        # operating conditions
        propemp.rpm = rpm
        propemp.ar_wing = ar_wing
        propemp.cl_wing = cl_wing
        propemp.cla_wing = cla_wing
        propemp.propulsor_type = propulsor_type
        propemp.power_condition = power_condition

        propemp.va_inlet = va_inlet
        propemp.d_exit = d_exit

        propemp.bem_input = bem_input

        propemp.delta_e = delta_e
        propemp.delta_r = delta_r

        propemp.conditions = conditions
        propemp.reference = reference
        propemp.density = air_density_isa(propemp.conditions[2])[0]
        propemp.comp_pe = comp_pe

        """ ------------------- determining the inflow angle based on wing downwash -------------------------------- """
        e0 = (2 * propemp.cl_wing) / (np.pi * propemp.ar_wing)
        de_da = 2 * propemp.cla_wing / (np.pi * propemp.ar_wing)
        eta = e0 + de_da * propemp.conditions[1]

        """ -------------------- conditions after the propeller disk ----------------------------------------------- """
        propemp.v_after_prop = np.sqrt((propemp.conditions[0] + 2 * propemp.bem_input[6]) ** 2
                                       + propemp.bem_input[7] ** 2)
        propemp.v_after_prop = 20
        propemp.a_after_prop = np.arctan(propemp.conditions[0] + 2 * propemp.bem_input[6] / propemp.bem_input[7])
        propemp.a_after_prop = 10

        # Initiate propeller class
        propemp.geometry_propeller = geometry_propeller
        propemp.propeller = Propeller(aircraft="DUUC",
                                      conditions=propemp.conditions, reference=propemp.reference,
                                      geometry=geometry_propeller, rpm=propemp.rpm,
                                      power_condition=propemp.power_condition, va_inlet=propemp.va_inlet,
                                      propulsor_type=propemp.propulsor_type, bem_input=propemp.bem_input,
                                      duct_diameter=geometry_duct[0])

        # Initiate duct class
        propemp.geometry_duct = geometry_duct
        propemp.duct = Duct(geometry=propemp.geometry_duct, conditions=propemp.conditions, reference=propemp.reference,
                            power_condition=propemp.power_condition, tc_prop=propemp.bem_input[3],
                            bem_input=propemp.bem_input, eta=eta)

        # Initiate nacelle class
        propemp.geometry_nacelle = geometry_nacelle
        propemp.nacelle = Nacelle(geometry=propemp.geometry_nacelle, conditions=propemp.conditions,
                                  reference=propemp.reference, propulsor_type=propemp.propulsor_type,
                                  power_condition=propemp.power_condition, v_after_prop=propemp.v_after_prop,
                                  a_after_prop=propemp.v_after_prop)

        # Initiate control vane class for elevator (1 piece)
        propemp.geometry_control = geometry_control
        propemp.elevator = ControlVane(geometry=propemp.geometry_control, conditions=propemp.conditions,
                                       reference=propemp.reference, power_condition=propemp.power_condition,
                                       v_after_prop=propemp.v_after_prop, a_after_prop=propemp.v_after_prop,
                                       va_inlet=propemp.va_inlet, d_exit=propemp.d_exit, deflection=propemp.delta_e)

        # Initiate control vane class for rudder (1 piece)
        propemp.rudder = ControlVane(geometry=propemp.geometry_control, conditions=propemp.conditions,
                                     reference=propemp.reference, power_condition=propemp.power_condition,
                                     v_after_prop=propemp.v_after_prop, a_after_prop=propemp.v_after_prop,
                                     va_inlet=propemp.va_inlet, d_exit=propemp.d_exit, deflection=propemp.delta_r)

        # Initiate support class
        propemp.geometry_support = geometry_support
        propemp.support = SupportStrut(geometry=propemp.geometry_support, conditions=propemp.conditions,
                                       reference=propemp.reference, power_condition=propemp.power_condition,
                                       prop_diameter=propemp.geometry_propeller[1], v_after_prop=propemp.v_after_prop,
                                       a_after_prop=propemp.v_after_prop, va_inlet=propemp.va_inlet)

        # Initiate pylon class
        propemp.geometry_pylon = geometry_pylon
        propemp.pylon = Pylon(geometry=propemp.geometry_pylon, conditions=propemp.conditions,
                              reference=propemp.reference, power_condition=propemp.power_condition)
    """ -------------------------------- geometric properties ------------------------------------------------------ """
    def inflow_angle(self):
        a_inflow = self.conditions[1]
        return a_inflow

    def area_wet(self):
        s_wet_duct = self.duct.area_wetted()
        s_wet_pylon = self.pylon.area_wetted()
        s_wet_support = self.support.area_wetted()
        s_wet_control = self.rudder.area_wetted()

        s_wet_total = 2 * (s_wet_duct + s_wet_pylon + s_wet_support) + 8 * s_wet_control
        return s_wet_total

    """ -------------------------------- coefficients -------------------------------------------------------------- """
    def cd0_vector(self):
        cd0_duct = self.duct.cd0()[1]
        cd0_pylon = self.pylon.cd0()
        cd0_nacelle = self.nacelle.cd0()[1]
        cd0_support = self.support.cd0()[1]
        cd0_control = self.elevator.cd0()[1] * 8
        return [cd0_duct, cd0_pylon, cd0_nacelle, cd0_control, cd0_control, cd0_support]

    def cd_interference_vector(self):
        cd_interference_pylon = self.pylon.cd_interference() * 2
        cd_interference_support = self.support.cd_interference() * 2
        cd_interference_control = self.rudder.cd_interference() * 8
        cd_interference_propeller = self.propeller.cd_interference() * 2

        return [cd_interference_pylon, cd_interference_support, cd_interference_control, cd_interference_propeller]

    def cd_vector(self):
        """ normalized vector"""
        cd_duct = self.duct.cd()[1]
        cd_pylon = self.pylon.cd()[1]
        cd_nacelle = self.nacelle.cd()[1]
        cd_support = self.support.cd()[1]
        cd_propeller = self.propeller.cd()[1]
        cd_control = 2 * (self.elevator.cd()[1] + self.rudder.cd()[1])
        return [cd_duct, cd_pylon, cd_nacelle, cd_support, cd_propeller, cd_control]

    def cl_sum(self):
        cl_sum_pe = (self.duct.cl()[1] + self.pylon.cl()[1] + self.support.cl()[1] + self.propeller.cl()[1]
                     + 2 * self.elevator.cl()[1])
        return cl_sum_pe

    def cd_sum(self):
        cd_sum_pe = (self.duct.cd()[1] + self.pylon.cd()[1] + self.support.cd()[1] + self.propeller.cd()[1]
                     + 2 * (self.elevator.cd()[1] + self.rudder.cd()[1]))
        return cd_sum_pe

    def cm_emp(self):
        """ assume the moment to be at the center line of the duct at the c/4 location """
        x_pylon = self.comp_pe[0]
        z_pylon = ((self.geometry_duct[0] / 2) * np.sin(np.radians(self.geometry_pylon[3]))
                   + (self.geometry_pylon[0] - self.geometry_duct[0] / 2) / 2
                   * np.sin(np.radians(self.geometry_pylon[3])))
        x_support = self.comp_pe[1]

        """ properties are all normalized with wing area, inflow velocity and cmac of the wing """
        cm_duct = self.duct.cm()[1]
        cm_support = self.support.cm()[1]
        cm_nac = self.nacelle.cm()[1]
        cm_pylon = self.pylon.cm()[1]

        support = x_support * self.support.cn()[1]
        pylon = x_pylon * self.pylon.cn()[1] + z_pylon * self.pylon.ct()[1]

        mf_sum = - support - pylon
        cm_sum = cm_duct + cm_nac + cm_pylon + cm_support

        propeller_moment = 0  # assume the propeller is located at c/4 and there is no moment in the x-z plane

        if self.power_condition == "off":
            cm_pe = cm_sum + mf_sum
            return cm_pe
        if self.power_condition == "on":
            cm_pe = cm_sum + mf_sum + propeller_moment
            return cm_pe
        else:
            raise ValueError("Empennage assembly PE.py -> wrong power condition input")

    def cl_norm_vector(self):
        cl_norm_duct = self.duct.cl()[1]
        cl_norm_pylon = self.pylon.cl()[1]
        cl_norm_nacelle = self.nacelle.cl()[1]
        cl_norm_support = self.support.cl()[1]
        cl_propeller = self.propeller.cl()[1]

        cl_norm_tot = cl_norm_duct + cl_norm_pylon + cl_norm_nacelle + cl_norm_support + cl_propeller
        return cl_norm_tot

    def cd_norm_vector(self):
        cd_norm_duct = self.duct.cd()[1]
        cd_norm_pylon = self.pylon.cd()[1]
        cd_norm_nacelle = self.nacelle.cd()[1]
        cd_norm_support = self.support.cd()[1]
        cd_propeller = self.propeller.cd()[1]

        cd_norm_tot = cd_norm_duct + cd_norm_pylon + cd_norm_nacelle + cd_norm_support + cd_propeller
        return cd_norm_tot

    """ -------------------------------------- Forces -------------------------------------------------------------- """
    def drag(self):
        drag_pe = (self.duct.drag_force() + self.pylon.drag_force() + self.support.drag_force()
                   + self.propeller.drag_force() + 2 * (self.rudder.drag_force() + self.elevator.drag_force()))
        return drag_pe

    def lift(self):
        lift_pe = (self.duct.lift_force() + self.pylon.lift_force() + self.support.lift_force()
                   + self.propeller.lift_force() + 2 * self.elevator.lift_force())
        return lift_pe

    def thrust(self):
        thrust_pe = 2 * (self.propeller.thrust() + self.duct.thrust())
        return thrust_pe

    def side_force(self):
        side_force = 0
        return side_force

    """ -------------------------------------- Weight ------------------------------------------------------------- """
    def weight_ps(self):
        # Safety factors
        n_ult = 1.5  # Ultimate load factor
        k_stoot = 1.5  # Impact factor
        k_misc = 1.25  # Miscellaneous factor

        # Convert cant angle to radians
        c_rad = np.radians(self.geometry_pylon[3])

        # Total length of the beam
        l_tot = self.geometry_support[0] + self.geometry_pylon[0]  # [m]

        """ Material properties -> AL 6061-T6 """
        sigma_allow = 241 * 1e6  # Allowable stress [Pa] (converted from MPa to N/m^2)
        rho = 2700  # Mass per unit length of the beam [density]

        # Aerodynamic forces on pylon and support
        f_pylon = (0.5 * self.density * self.pylon.inflow_velocity() ** 2 * self.pylon.cl()[0] * self.geometry_pylon[1]
                   * self.geometry_pylon[0])  # [N]
        f_support = (0.5 * self.density * self.support.inflow_velocity() ** 2 * self.support.cl()[0]
                     * self.geometry_support[1] * self.geometry_support[0])  # [N]

        # Aerodynamic moment about the root
        m_aero = (0.5 * f_pylon * self.geometry_pylon[0] + (0.5 * self.geometry_support[0] + self.geometry_pylon[0]) * f_support) * np.cos(
            c_rad)  # [N·m]

        # Weight moment about the root
        m_weight = (self.duct.weight() + self.nacelle.weight() + self.propeller.weight_engine()
                    + self.propeller.weight_fan()) * (self.geometry_pylon[0] + self.geometry_support[0] * 0.5) * 9.81 * np.cos(
            c_rad)  # [N·m]

        # Sizing moment (maximum of aerodynamic and weight moments)
        m_sizing = max(m_aero, m_weight)  # [N·m]

        # Determine maximum thickness of the airfoil section
        num_list = [int(digit) for digit in self.geometry_pylon[2]]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness percentage
        thickness_pylon = (thickness / 100) * self.geometry_pylon[1]  # [m]

        num_list2 = [int(digit) for digit in self.geometry_support[2]]
        thickness_sup = num_list2[2] * 10 + num_list2[3]  # NACA thickness percentage
        thickness_support = (thickness_sup / 100) * self.geometry_support[1]  # [m]

        height = max(thickness_pylon, thickness_support) * 0.90  # Reduce maximum thickness by a safety factor [m]

        # Calculate width using bending stress formula
        width = (((6 * m_sizing * n_ult * k_stoot) / sigma_allow) / (height ** 2)) * k_misc  # [m]

        # Total mass of the beam
        mass_tot = height * width * l_tot * rho  # [kg]

        # Split mass into pylon and support sections
        m_pylon = mass_of_section(mass_tot, l_tot, 0, self.geometry_pylon[0])  # Mass of pylon section [kg]
        m_support = mass_of_section(mass_tot, l_tot, self.geometry_pylon[0], l_tot)  # Mass of support section [kg]

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
        alpha = np.radians(self.conditions[1])

        cn_dp_nasa = f1 * np.sin(alpha) * (np.cos(alpha) + f2 * self.yv())
        return cn_dp_nasa

    def ct_dp(self):
        f3 = 1.90
        f4 = .92
        alpha = np.radians(self.conditions[1])

        ct_dp_nasa = f3 * np.sin(alpha) ** 2 + f4 * (self.yv()) ** 2
        return ct_dp_nasa

    def yv(self):
        alpha = np.radians(self.conditions[1])
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
        alpha = np.radians(self.conditions[1])
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

    cl_norm = []
    cd_norm = []

    for i in range(len(a)):
        PE = PropulsiveEmpennage(rpm=config.rpm,
                                 alpha=a[i],
                                 power_condition="off",
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
                                 bem_input=[4000., 26482, -1.4475057750305072e-13, 0.8892, 0.3297, 0.3201, 5, 10, 1820],
                                 delta_e=0, delta_r=0, altitude=7000)
        cd_intereference.append(PE.cd_interference_vector())
        cl.append(PE.cl_or())
        cd.append(PE.cd_or())
        cd2.append(PE.cd_sum())
        cm.append(PE.cm_emp())
        cl_norm.append(PE.cl_norm())
        cd_norm.append(PE.cd_norm())

        polar_pylon = airfoil_polar(f"pylon0012.txt", float(a[i]))
        cd_val_pylon = float(polar_pylon[1])
        cl_val_pylon = float(polar_pylon[0])
        cm_val_pylon = float(polar_pylon[2])

        polar_support = airfoil_polar(f"support0012.txt", float(a[i]))
        cd_val_support = float(polar_support[1])
        cl_val_support = float(polar_support[0])
        cm_val_support = float(polar_support[2])

        ct_dp.append(PE.ct_dp())
        cn_dp.append(PE.cn_dp())
        cm_dp.append(PE.cm_dp())



    alpha_avl, cl_avl, clff_avl, cd_avl, cdin_avl, cdff_avl = read_avl_output("AVL_RW_aeroproperties.txt")

    results = np.array(cd_intereference)
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
    """

    plt.figure('Cl - alpha')
    plt.plot(a, cl, label=r'Prediction model', color="tab:blue")
    plt.plot(refavl.a, [cl / (10 / 1.5) for cl in refavl.cl_pe_nopylon], label=r'AVL', color='tab:purple', linestyle='dashed', marker='x')
    plt.plot(ref7f.pr_cla_a, ref7f.pr_cla_cl, label='Experimental', color='tab:orange', linestyle='dashed')
    plt.plot(a, cn_dp, label='Nasa prediction', color='red', linestyle='dashed')
    plt.xlim([0, 20])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{l}$ [-]')
    plt.title(r'$C_{l}$ vs. $\alpha$ - Propulsive Empennage')
    plt.legend()
    plt.grid(True)

    plt.figure('Cd - alpha')
    plt.plot(a, cd, label=r'Prediction model', color="tab:blue")
    plt.plot(refavl.a, [0.0175 +cd / (10 / 1.5) for cd in refavl.cd_pe_nopylon], label=r'AVL', color='tab:purple', linestyle='dashed', marker='x')
    #plt.plot(ref7f.pr_cda_a, ref7f.pr_cla_cl, label='Experimental', color='tab:orange', linestyle='dashed')
    plt.plot(a, ct_dp, label='Nasa prediction', color='red', linestyle='dashed')
    plt.xlim([0, 20])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{d}$ [-]')
    plt.title(r'$C_{d}$ vs. $\alpha$ - Propulsive Empennage')
    plt.legend()
    plt.grid(True)


    plt.figure('CL - CD')

    # Plot Prediction model
    plt.plot(cd, cl, label=r'Prediction model', color="tab:blue")

    # Get AVL data
    avl_cd = [0.0176 + (cd / (10 / 1.5)) for cd in refavl.cd_pe_nopylon]
    avl_cl = [cl / (10 / 1.5) for cl in refavl.cl_pe_nopylon]

    # Get Experimental data
    exp_cd = [cd - 0.07 for cd in ref7f.pr_clcd_cd]
    exp_cl = ref7f.pr_clcd_cl

    # Sort both datasets by Cd for proper interpolation/extrapolation
    avl_sorted = sorted(zip(avl_cd, avl_cl))
    exp_sorted = sorted(zip(exp_cd, exp_cl))

    avl_cd_sorted, avl_cl_sorted = zip(*avl_sorted)
    exp_cd_sorted, exp_cl_sorted = zip(*exp_sorted)

    # Create interpolation functions with extrapolation enabled
    avl_interp_func = interp1d(avl_cd_sorted, avl_cl_sorted, kind='linear', fill_value="extrapolate")
    exp_interp_func = interp1d(exp_cd_sorted, exp_cl_sorted, kind='linear', fill_value="extrapolate")

    # Define extended Cd range for shading (slightly beyond bounds)
    extended_cd_min = min(min(avl_cd_sorted), min(exp_cd_sorted)) - 0.01
    extended_cd_max = max(max(avl_cd_sorted), max(exp_cd_sorted)) + 0.01
    common_cd = np.linspace(extended_cd_min, extended_cd_max, 500)

    # Interpolate and extrapolate Cl values for both curves
    avl_cl_interp = avl_interp_func(common_cd)
    exp_cl_interp = exp_interp_func(common_cd)

    # Plot original curves
    plt.plot(avl_cd_sorted, avl_cl_sorted, label=r'AVL', color='tab:purple', marker='x')
    plt.plot(exp_cd_sorted, exp_cl_sorted, label='Experimental', color='tab:orange', linestyle='dashed')

    # Add shaded area between interpolated/extrapolated curves
    plt.fill_between(common_cd, avl_cl_interp, exp_cl_interp, color='tab:blue', alpha=0.5)

    # Add labels and formatting
    plt.xlabel(r'$C_{d}$ [-]')
    plt.ylabel(r'$C_{l}$ [-]')
    plt.xlim([0, 0.25])
    plt.ylim([0, 2])
    plt.title(r'$C_{l}$ vs. $C_{d}$ - Propulsive Empennage')
    plt.legend(loc='lower right')
    plt.grid(True)

    plt.figure('CL - CD (normalized)')
    plt.plot(cd_norm, cl_norm, label=r'Prediction model', color="tab:blue")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $C_{D}$ - Propulsive Empennage')
    plt.legend()
    plt.grid(True)
    plt.show()

"""
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
    plt.show()
"""

