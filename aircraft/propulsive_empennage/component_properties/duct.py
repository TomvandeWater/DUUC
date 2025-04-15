from analysis_modules.aerodynamic import reynolds
from analysis_modules.ISA import air_density_isa
from analysis_modules.factors import *
from data.read_data import *
import data.atr_reference as ref
import config
import matplotlib.pyplot as plt
import data.experiment_reference_7foot as ref7f
import data.experiment_reference_5annular_airfoil as ref5r
import data.experiment_avl as refavl
from scipy.interpolate import interp1d
import data.experiment_reference_linear_regression as reflr


class Duct:
    def __init__(self, geometry, conditions, reference, power_condition: str, tc_prop: float, bem_input,
                 eta: float):
        super().__init__()
        self.duct_diameter = geometry[0]
        self.duct_chord = geometry[1]
        self.duct_profile = geometry[2]

        self.v_inf = conditions[0]
        self.alpha = conditions[1]
        self.altitude = conditions[2]
        self.mach = conditions[3]
        self.density = air_density_isa(self.altitude)[0]
        self.a_install_wing = conditions[4]
        self.a_install_duct = conditions[5]
        self.beta = conditions[6]

        self.pc = power_condition

        self.ref_area = reference[0]
        self.ref_chord = reference[1]

        self.bem_input = bem_input
        self.va = self.bem_input[8]
        self.tc_prop = tc_prop
        self.eta = eta

    """ --------------------------------------- Define inflow properties --------------------------------------- """
    def inflow_velocity(self):
        if self.pc == "off":
            u_duct = self.v_inf
            return u_duct
        else:
            u_duct = self.v_inf
            return u_duct

    def inflow_angle(self):
        """ inflow angle for the duct is equal to the freestream velocity"""
        inflow_duct = self.alpha - self.a_install_wing - self.eta + self.a_install_duct
        angle_rad = np.deg2rad(inflow_duct)
        return inflow_duct, angle_rad

    def reynolds_number(self):
        """ local reynolds number of the duct, velocity is equal to the velocity inside the duct"""
        re_duct = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.duct_chord)
        return re_duct

    """ ----------------------------------------- Determine geometric properties --------------------------------- """
    def area_wetted(self):
        num_list = [int(digit) for digit in self.duct_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        t_max = thickness / 100

        s_wet = (2 * np.pi * self.duct_chord * self.duct_diameter + 2 * np.pi
                 * self.duct_diameter * 0.5 * t_max * self.duct_chord)
        return s_wet

    def proj_area(self):
        """ projected area of the duct"""
        proj_area = self.duct_chord * self.duct_diameter
        return proj_area

    def aspect_ratio(self):
        """ aspect ratio defined as diameter divided by chord -> note deviation with standard wings"""
        ar_duct = self.duct_diameter / self.duct_chord
        return ar_duct

    def t_c(self):
        """ assume NACA 44-series airfoil"""
        num_list = [int(digit) for digit in self.duct_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    """ ------------------------------ Ratios used for normalization --------------------------------------------- """
    def area_ratio(self):
        ar_pylon = self.proj_area() / self.ref_area
        return ar_pylon

    def area_ratio_wet(self):
        ar_w_pylon = self.area_wetted() / self.ref_area
        return ar_w_pylon

    def velocity_ratio(self):
        v_ratio = self.inflow_velocity() ** 2 / self.v_inf ** 2
        return v_ratio

    def chord_ratio(self):
        c_ratio = self.duct_chord / self.ref_chord
        return c_ratio

    """" -------------------------------------- coefficient ------------------------------------------------------ """
    def zeta(self):
        """ based on Weissinger prediction model"""
        delta = 1 / self.aspect_ratio()
        zeta_duct = 1 / (1 + delta * np.pi / 2 + np.arctan(1.2 * delta) * delta)
        return zeta_duct

    @staticmethod
    def cl_a():
        cl_a_naca = 2 * np.pi
        return cl_a_naca

    def cl_da(self):
        cl_da_duct = np.pi / 2 * self.zeta() * self.cl_a()
        return cl_da_duct

    def cl(self):
        k_prop = 0.2 * np.sqrt(np.abs(self.tc_prop))
        norm_area = self.area_ratio()
        norm_speed = self.velocity_ratio()

        if self.pc == "off":
            cl_duct = self.cl_da() * self.inflow_angle()[1]

            cl_duct_norm = cl_duct * norm_area * norm_speed

        elif self.pc == "on":
            cl_duct = (1 + k_prop) * self.cl_da() * self.inflow_angle()[1]

            cl_duct_norm = cl_duct * norm_area * norm_speed
        else:
            raise ValueError("Power conditions not specified properly")

        return cl_duct, cl_duct_norm

    def cy(self):
        k_prop = 0.2 * np.sqrt(np.abs(self.tc_prop))
        norm_area = self.area_ratio()
        norm_speed = self.velocity_ratio()

        if self.pc == "off":
            cy_duct = self.cl_da() * self.beta
            cy_duct_norm = cy_duct * norm_area * norm_speed

        elif self.pc == "on":
            cy_duct = (1 + k_prop) * self.cl_da() * self.beta
            cy_duct_norm = cy_duct * norm_area * norm_speed
        else:
            raise ValueError("Power conditions not specified properly")

        return cy_duct, cy_duct_norm


    def cd0(self):
        cf = skin_friction(self.reynolds_number(), 't')
        fm = mach_correction(self.mach)
        ftc = 1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4

        coeff = airfoil_polar(f"duct{self.duct_profile}.txt", float(0.0))
        cdmin = float(coeff[1])

        cd0_norm = fm * ftc * cf * self.area_ratio_wet() #* (cdmin / 0.004) ** 4
        cd0_duct = fm * ftc * cf * (cdmin / 0.004) ** 4
        return cd0_duct, cd0_norm

    def cdi(self):
        oswald_duct = 2  # from vikesh

        cdi_duct = (self.cl()[0]) ** 2 / (oswald_duct * np.pi * self.aspect_ratio())

        cdi_norm = (((self.cl()[0] ** 2) / (oswald_duct * np.pi * self.aspect_ratio())) * self.area_ratio()
                    * self.velocity_ratio())
        return cdi_duct, cdi_norm

    def cd(self):
        cd_duct = self.cdi()[0] + self.cd0()[0]
        cd_duct_norm = self.cdi()[0] + self.cd0()[1] * self.velocity_ratio()
        return cd_duct, cd_duct_norm

    def cm(self):
        """ Aerodynamic moment of an annular wing based on Masqood"""
        alfa = self.inflow_angle()[1]
        ac = 0.3
        b = -0.786
        c = - 0.26

        p1 = 0.63
        p0 = -0.35

        kpcm = 6.25 * np.sin(self.aspect_ratio()/2)
        kvcm = np.pi / 3

        xp = ac * self.aspect_ratio() ** b + c
        xe = p1 * self.aspect_ratio() + p0

        cm_duct = xp * kpcm * np.sin(alfa) * np.cos(alfa) + xe * kvcm * np.sin(alfa) ** 2
        cm_duct_norm = cm_duct * self.area_ratio() * self.velocity_ratio() * self.chord_ratio()

        return cm_duct, cm_duct_norm

    def cn(self):
        alpha = self.inflow_angle()[1]

        cn_duct = self.cl()[0] * np.cos(alpha) + self.cd()[0] * np.sin(alpha)

        cn_norm = cn_duct * self.area_ratio() * self.velocity_ratio()
        return cn_duct, cn_norm

    def ct(self):
        alpha = self.inflow_angle()[1]

        ct_duct = self.cl()[0] * np.sin(alpha) - self.cd()[0] * np.cos(alpha)

        ct_norm = ct_duct * self.area_ratio() * self.velocity_ratio()
        return ct_duct, ct_norm

    """ --------------------------------------------------- FORCES ---------------------------------------------- """
    def thrust(self):
        """ This thrust component is only dependent on the shape of the duct, the pressure induced thrust on the
        leading edge is not incorporated here"""
        area_prop = np.pi * (ref.blade_diameter / 2) ** 2 - (np.pi * (config.nacelle_diameter / 2) ** 2)
        area_exit = np.pi * (config.d_exit / 2) ** 2 - (np.pi * (config.nacelle_diameter / 2) ** 2)

        v_exit = ((self.v_inf + 2 * self.va) * area_prop) / area_exit
        m_dot = self.density * area_prop * (v_exit + self.v_inf) / 2

        thrust_tot = m_dot * (v_exit - self.v_inf)
        t_duct = thrust_tot - self.bem_input[0]

        ct_duct = t_duct / (0.5 * self.density * self.v_inf ** 2 * self.ref_area)

        return t_duct, ct_duct

    def lift_force(self):
        lift_duct = self.cl()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.proj_area()
        return lift_duct

    def drag_force(self):
        drag_duct = self.cd()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.proj_area()
        return drag_duct

    def moment_force(self):
        moment_duct = (self.cm()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.proj_area()
                       * self.duct_chord)
        return moment_duct

    """" ------------------------------------ NASA PREDICTION MODEL ------------------------------------------------ """
    def cn_dp(self):
        f1 = 3.10
        f2 = .53
        alpha = np.radians(self.inflow_angle())

        cn_dp_nasa = f1 * np.sin(alpha) * (np.cos(alpha) + f2 * self.yv())
        return cn_dp_nasa

    def ct_dp(self):
        f3 = 1.90
        f4 = .92
        alpha = np.radians(self.inflow_angle())

        ct_dp_nasa = f3 * np.sin(alpha) ** 2 + f4 * (self.yv()) ** 2
        return ct_dp_nasa

    def yv(self):
        alpha = np.radians(self.inflow_angle())
        f4 = .92
        f3 = 1.90
        C_T_dp_guess = 0.168

        A = (config.d_exit / 2) ** 2 * np.pi
        A_p = (config.d_exit / 2) ** 2 * np.pi - (config.nacelle_diameter / 2) ** 2 * np.pi

        a_r = A / A_p

        a_c = -1 * np.cos(alpha)/(a_r * f4 + 1)
        b = np.cos(alpha)/(a_r * f4 + 1)
        d = a_r * C_T_dp_guess - (a_r * f3 - 1) * (np.sin(alpha) ** 2)
        c = a_r * f4 + 1

        yv_cal = a_c + np.sqrt(b**2 + d / c)
        return yv_cal

    def cm_dp(self):
        alpha = np.radians(self.inflow_angle())
        f5 = .22
        f6 = 1.5
        f7 = 1.87

        yv = self.yv()

        cm_pe = 4 * f5 * np.sin(alpha) * np.cos(alpha) + (f5 * f6 + f7) * yv * np.sin(alpha)

        return cm_pe

    """ -------------------------------------------- Weights --------------------------------------------------- """
    def weight(self):
        """ based on Torenbeek class II weight estimation for horizontal stabilizer """
        kh = 1.05
        sh = self.duct_diameter * np.pi * self.duct_chord
        vd = ref.v_dive
        sweep = np.radians(0)

        m_duct = kh * sh * (62 * (sh ** 0.2 * vd) / (1000 * np.sqrt(np.cos(sweep))) - 2.5)
        return m_duct


""" Test section"""
"""
if __name__ == "__main__":

    a = np.linspace(0, 20, 41)
    cl = []
    cd = []
    cm = []
    cl_the = []
    cd_the = []
    for i in range(len(a)):
        conditions = [128, a[i], 7000, 0.41, 0, 0]
        geometry = [config.duct_diameter, config.duct_chord, "0012"]
        reference = [62, 2.5]
        wing = Duct(conditions=conditions, geometry=geometry,
                    power_condition="off",
                    tc_prop=0.48,
                    bem_input=[5000, 26482.06279555917, 1800, -1.4475057750305072e-13, 0.8892292886261024, 0.3297344029147765, 0.3201225439053968, 5, 10],
                    eta=0, reference=reference)
        al = np.radians(a[i])
        kp = 6.25 * np.sin(wing.aspect_ratio()/2)
        kv = np.pi / 3
        cl_theory = kp * np.sin(al) * np.cos(al)**2 + kv * np.cos(al) * np.sin(al) ** 2
        cl_the.append(cl_theory)
        cd_the.append(0.0125 + 0.06 * cl_theory ** 2)

        cl.append(wing.cl()[0])
        cd.append(wing.cd()[0])
        cm.append(wing.cm()[0])

    alpha_avl, cl_avl, clff_avl, cd_avl, cdin_avl, cdff_avl = read_avl_output("AVL_RW_aeroproperties.txt")

    plt.figure('CL - alpha')
    plt.plot(a, cl, label=r'Prediction model', color="tab:blue")
    plt.plot(ref5r.cla_a_ar_1_5, ref5r.cla_cl_ar_1_5, label=r'Experimental', color='tab:green', linestyle='dashed', marker='o')
    plt.plot(refavl.a, [cl_rw / (10 / wing.aspect_ratio()) for cl_rw in refavl.cl_rw], label=r'AVL', color='tab:purple', linestyle='dashed', marker='x')
    plt.plot(a, cl_the, label='Leading Edge Suction analogy', color="tab:red", linestyle="--")
    plt.xlim([0, 20])
    plt.ylim([-0.1, 1.7])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{l}$ [-]')
    plt.title(r'$C_{l}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    cl_the_interp = np.interp(a, ref5r.cla_a_ar_1_5, ref5r.cla_cl_ar_1_5)

    plt.figure('CL - alpha (version 2)')
    plt.plot(a, cl, label=r'Prediction model', color="tab:blue")
    plt.fill_between(a, cl_the, cl_the_interp, where=(cl_the >= cl_the_interp), interpolate=True, color="tab:blue", alpha=0.5,)
    plt.plot(ref5r.cla_a_ar_1_5, ref5r.cla_cl_ar_1_5, label=r'Experimental', color='tab:green', linestyle='dashed')
    plt.plot(a, cl_the, label='Leading Edge Suction analogy', color="tab:red", linestyle="--")
    plt.xlim([0, 20])
    plt.ylim([-0.1, 1.7])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{l}$ [-]')
    plt.title(r'$C_{l}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - alpha (version 1)')
    plt.plot(a, cl, label=r'Prediction model', color="tab:blue")
    plt.fill_between(a, cl_the, cl_the_interp, where=(cl_the >= cl_the_interp), interpolate=True, color="tab:blue", alpha=0.5)
    plt.xlim([0, 20])
    plt.ylim([-0.1, 1.7])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{l}$ [-]')
    plt.title(r'$C_{l}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - alpha (normalized)')
    plt.plot(a, [cl * (wing.proj_area() / wing.ref_area) for cl in cl], label=r'Prediction model', color="tab:blue")
    plt.xlim([0, 20])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Prediction model', color="tab:blue")
    plt.plot(ref5r.cda_a_ar_1_5, ref5r.cda_cd_ar_1_5, label=r'Experimental', color='tab:green', linestyle='dashed', marker='o')
    plt.plot(refavl.a, [0.0125 + cd_rw / (10 / wing.aspect_ratio()) for cd_rw in refavl.cd_rw], label=r'AVL', color='tab:purple', linestyle='dashed', marker='x')
    plt.plot(a, cd_the, label='Leading Edge Suction analogy', color="tab:red", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $\alpha$ - Duct')
    plt.xlim([0, 20])
    plt.ylim([0, 0.25])
    plt.legend()
    plt.grid(True)

    cd_the_interp = np.interp(a, ref5r.cda_a_ar_1_5, ref5r.cda_cd_ar_1_5)

    plt.figure('CD - alpha (version 2)')
    plt.plot(a, cd, label=r'Prediction model', color="tab:blue")
    plt.fill_between(a, cd_the_interp, cd_the, where=(cd_the_interp >= cd_the), interpolate=True, color="tab:blue", alpha=0.5,)
    plt.plot(ref5r.cda_a_ar_1_5, ref5r.cda_cd_ar_1_5, label=r'Experimental', color='tab:green', linestyle='dashed')
    plt.plot(a, cd_the, label='Leading Edge Suction analogy', color="tab:red", linestyle="--")
    plt.xlim([0, 20])
    plt.ylim([0, 0.25])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{d}$ [-]')
    plt.title(r'$C_{d}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('Cd - alpha (version 1)')
    plt.plot(a, cd, label=r'Prediction model', color="tab:blue")
    plt.fill_between(a, cd_the_interp, cd_the, where=(cd_the_interp >= cd_the), interpolate=True, color="tab:blue", alpha=0.5,)
    plt.xlim([0, 20])
    plt.ylim([0, 0.25])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{d}$ [-]')
    plt.title(r'$C_{d}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha (normalized)')
    plt.plot(a, [cd * (wing.proj_area() / wing.ref_area) for cd in cd], label=r'Prediction model', color="tab:blue")
    plt.xlim([0, 20])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CM - alpha')
    plt.plot(a, cm, label=r'Prediction model', color="tab:blue")
    #plt.plot(a_ref, cm_ref, label=r'Experimental', color="tab:green", marker='o')
    plt.plot(ref5r.cma_a_ar_1_5, ref5r.cma_cm_ar_1_5, label=r'Experimental', color='tab:green', linestyle='dashed')
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{M}$ [-]')
    plt.title(r'$C_{M}$ vs. $\alpha$ - Duct')
    plt.xlim([0, 20])
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD - validation')
    plt.plot(cd, cl, label=r'Prediction model', color="tab:blue")
    plt.plot([0.0125 + cla ** 2 / (2 * np.pi * 1.5) for cla in ref5r.cla_cl_ar_1_52], ref5r.cla_cl_ar_1_52, label=r'Experimental', color='tab:green', linestyle='dashed')
    plt.plot(cd_the, cl_the, label='Leading Edge Suction analogy', color="tab:red", linestyle="--")
    plt.xlabel(r'$C_{d}$ [-]')
    plt.ylabel(r'$C_{l}$ [-]')
    plt.xlim([0, 0.25])
    plt.ylim([-0.1, 1.6])
    plt.title(r'$C_{l}$ vs. $C_{d}$ - Duct')
    plt.legend()
    plt.grid(True)

    cd_exp = [0.0125 + cla ** 2 / (2 * np.pi * 1.5) for cla in ref5r.cla_cl_ar_1_52]
    cl_exp = ref5r.cla_cl_ar_1_52

    # Interpolation functions for both lines
    interp_exp = interp1d(cl_exp, cd_exp, kind='linear', fill_value="extrapolate")
    interp_the = interp1d(cl_the, cd_the, kind='linear', fill_value="extrapolate")

    # Define the range of C_L values for shading (only for C_L > 0.05)
    cl_shade = np.linspace(0.05, min(max(cl_exp), max(cl_the)), 500)

    # Get corresponding C_D values for both lines
    cd_shade_exp = interp_exp(cl_shade)
    cd_shade_the = interp_the(cl_shade)

    # Plotting
    plt.figure('CL - CD - (version 1)')
    plt.plot(cd, cl, label=r'Prediction model', color="tab:blue")
    plt.plot(cd_exp, cl_exp, label=r'Experimental', color='tab:green', linestyle='dashed')
    plt.plot(cd_the, cl_the, label='Leading Edge Suction analogy', color="tab:red", linestyle="--")

    # Shade the area between the two lines
    plt.fill_betweenx(cl_shade, cd_shade_exp, cd_shade_the, color='tab:blue', alpha=0.5)
    plt.xlabel(r'$C_{d}$ [-]')
    plt.ylabel(r'$C_{l}$ [-]')
    plt.xlim([0, 0.25])
    plt.ylim([-0.1, 1.6])
    plt.title(r'$C_{l}$ vs. $C_{d}$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD - (version 2)')
    plt.plot(cd, cl, label=r'Prediction model', color="tab:blue")
    plt.fill_betweenx(cl_shade, cd_shade_exp, cd_shade_the, color='tab:blue', alpha=0.5)
    plt.xlabel(r'$C_{d}$ [-]')
    plt.ylabel(r'$C_{l}$ [-]')
    plt.xlim([0, 0.25])
    plt.ylim([-0.1, 1.6])
    plt.title(r'$C_{l}$ vs. $C_{d}$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD - (normalized)')
    plt.plot([cd * (wing.proj_area() / wing.ref_area) for cd in cd], [cl * (wing.proj_area() / wing.ref_area) for cl in cl], label=r'Prediction model', color="tab:blue")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $C_{D}$ - Duct')
    plt.legend()
    plt.grid(True)

    cl_conv = [4.999999999999449e-05, 0.03589707460042309, 0.07174414920084618, 0.10759122380126929, 0.14343829840169237,
     0.17928537300211547, 0.21513244760253858, 0.25097952220296166, 0.28682659680338474, 0.32267367140380787,
     0.35852074600423095, 0.39436782060465403, 0.43021489520507716, 0.46606196980550024, 0.5019090444059233,
     0.5377561190063463, 0.5736031936067695, 0.6094502682071927, 0.6452973428076157, 0.6811444174080389,
     0.7169914920084619, 0.752838566608885, 0.7886856412093081, 0.8245327158097312, 0.8603797904101543,
     0.8962268650105774, 0.9320739396110005, 0.9679210142114235, 1.0037680888118468, 1.03961516341227,
     1.0754622380126926, 1.1113093126131157, 1.1471563872135389, 1.183003461813962, 1.2188505364143856,
     1.2546976110148083, 1.2905446856152314, 1.3263917602156545, 1.3622388348160777, 1.3980859094165008,
     1.433932984016924]
    cd_conv =[0.008492249890895883, 0.008646092674148053, 0.009106765083414481, 0.009874267118695169, 0.010948598779990116,
     0.012329760067299323, 0.014017750980622789, 0.016012571519960512, 0.018314221685312493, 0.02092270147667874,
     0.023838010894059243, 0.027060149937454002, 0.030589118606863024, 0.0344249169022863, 0.038567544823723834,
     0.043017002371175624, 0.04777328954464169, 0.052836406344122025, 0.05820635276961659, 0.06388312882112543,
     0.0698667344986485, 0.07615716980218587, 0.08275443473173746, 0.08965852928730333, 0.09686945346888348,
     0.10438720727647785, 0.1122117907100865, 0.12034320376970939, 0.1287814464553466, 0.13752651876699803,
     0.14657842070466362, 0.15593715226834357, 0.1656027134580378, 0.17557510427374626, 0.18585432471546912,
     0.19644037478320597, 0.20733325447695725, 0.21853296379672274, 0.2300395027425025, 0.24185287131429653,
     0.25397306951210485]
    cl_norm_conv = [9.534276437384782e-06, 0.0068450526506779056, 0.013680571024918426, 0.02051608939915895, 0.02735160777339947,
     0.034187126147639996, 0.04102264452188052, 0.047858162896121034, 0.05469368127036156, 0.06152919964460209,
     0.0683647180188426, 0.07520023639308313, 0.08203575476732365, 0.08887127314156416, 0.09570679151580469,
     0.10254230989004519, 0.10937782826428573, 0.11621334663852627, 0.12304886501276678, 0.1298843833870073,
     0.13671990176124782, 0.14355542013548833, 0.15039093850972887, 0.15722645688396938, 0.1640619752582099,
     0.17089749363245044, 0.17773301200669095, 0.18456853038093146, 0.191404048755172, 0.19823956712941251,
     0.20507508550365297, 0.2119106038778935, 0.21874612225213405, 0.2255816406263746, 0.23241715900061516,
     0.23925267737485562, 0.24608819574909616, 0.2529237141233367, 0.2597592324975772, 0.2665947508718178,
     0.27343026924605823]
    cd_norm_conv = [0.0028944859851903234, 0.0029238215776587837, 0.0030116651395992373, 0.003158016671011684, 0.003362876171896124,
     0.003626243642252557, 0.0039481190820809834, 0.004328502491381402, 0.0047673938701538136, 0.005264793218398221,
     0.0058207005361146195, 0.006435115823303011, 0.007108039079963397, 0.007839470306095775, 0.008629409501700144,
     0.009477856666776509, 0.010384811801324867, 0.011350274905345221, 0.012374245978837564, 0.013456725021801902,
     0.014597712034238229, 0.015797207016146553, 0.01705520996752687, 0.01837172088837918, 0.019746739778703486,
     0.02118026663849978, 0.02267230146776807, 0.02422284426650835, 0.02583189503472063, 0.027499453772404897,
     0.029225520479561142, 0.0310100951561894, 0.03285317780228965, 0.03475476841786189, 0.03671486700290615,
     0.038733473557422356, 0.04081058808141059, 0.042946210574870795, 0.04514034103780301, 0.04739297947020721,
     0.049704125872083414]

    plt.figure('Duct vs. Empennage')
    plt.plot(cd, cl, label=r'Duct', color="tab:blue")
    plt.plot(cd_conv, cl_conv, label='Conventional Empennage', color='tab:orange')
    plt.xlabel(r'$C_{d}$ [-]')
    plt.ylabel(r'$C_{l}$ [-]')
    plt.xlim([0, 0.25])
    plt.ylim([-0.1, 1.6])
    plt.title(r'$C_{l}$ vs. $C_{d}$ - Comparison')
    plt.legend()
    plt.grid(True)

    plt.figure('Duct vs. Empennage (normalized)')
    plt.plot([cd * (wing.proj_area() / wing.ref_area) for cd in cd], [cl * (wing.proj_area() / wing.ref_area) for cl in cl], label=r'Duct', color="tab:blue")
    plt.plot(cd_norm_conv, cl_norm_conv, label='Conventional Empennage', color='tab:orange')
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $C_{D}$ - Comparison')
    plt.legend()
    plt.grid(True)

    plt.show()

    plt.figure('CD - CL^2')
    plt.plot([cl_val ** 2 for cl_val in cl], cd, label=r'Prediction model', color="tab:blue")
    #plt.plot([cl_val ** 2 for cl_val in cl_the], cd_the, label='Theory', color="tab:green", linestyle="--")
    #plt.plot([cl_val ** 2 for cl_val in cl], cd_0, label=r'$C_{D0}$', color="tab:blue", linestyle="--")
    #plt.plot([cl_avl ** 2 for cl_avl in cl_avl], cd_avl, label="AVL", color="purple")
    plt.xlabel(r'$C_{L}^{2}$ [-]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $C_{L}^{2}$ - Duct')
    plt.xlim([0, 0.5])
    plt.ylim([0, 0.04])
    plt.legend()
    plt.grid(True)
    plt.show()
    """
