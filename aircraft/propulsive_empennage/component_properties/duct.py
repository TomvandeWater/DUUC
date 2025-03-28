from analysis_modules.aerodynamic import reynolds
from analysis_modules.ISA import air_density_isa
from analysis_modules.factors import *
from data.read_data import *
import data.atr_reference as ref
import config
import matplotlib.pyplot as plt


class Duct:
    def __init__(self, duct_diameter: float, duct_chord: float, duct_profile: str,
                 alpha: float, power_condition: str,
                 tc_prop: float, v_inf: float, mach: float, ref_area: float, ref_chord: float, bem_input, va: float,
                 altitude: float):
        super().__init__()
        self.duct_diameter = duct_diameter
        self.duct_chord = duct_chord
        self.duct_profile = duct_profile
        self.alpha = alpha
        self.pc = power_condition
        self.tc_prop = tc_prop
        self.v_inf = v_inf
        self.mach = mach
        self.ref_area = ref_area
        self.ref_chord = ref_chord
        self.bem_input = bem_input
        self.va = va
        self.altitude = altitude
        self.density = air_density_isa(self.altitude)[0]
    """ --------------------------------------- Define inflow properties --------------------------------------- """
    def inflow_velocity(self):
        if self.pc == "off":
            u_duct = self.v_inf
            return u_duct
        else:
            u_duct = self.v_inf
            return u_duct

    def inflow_angle(self):
        inflow_duct = self.alpha
        return inflow_duct

    def reynolds_number(self):
        re_duct = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.duct_chord)
        return re_duct

    """ ----------------------------------------- Determine geometric properties --------------------------------- """
    def wetted_area(self):
        num_list = [int(digit) for digit in self.duct_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        t_max = thickness / 100

        s_wet = (2 * np.pi * self.duct_chord * self.duct_diameter + 2 * np.pi
                 * self.duct_diameter * 0.5 * t_max * self.duct_chord)
        return s_wet

    def proj_area(self):
        proj_area = self.duct_chord * self.duct_diameter
        return proj_area

    def aspect_ratio(self):
        ar_duct = self.duct_diameter / self.duct_chord
        return ar_duct

    def t_c(self):
        num_list = [int(digit) for digit in self.duct_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    """" -------------------------------------- coefficient ------------------------------------------------------ """
    def zeta(self):
        """ based on Weissinger prediction model"""
        delta = 1 / self.aspect_ratio()
        zeta_duct = 1 / (1 + delta * np.pi / 2 + np.arctan(1.2 * delta) * delta)
        return zeta_duct

    def cl_a(self):
        cl_a_naca = 2 * np.pi
        return cl_a_naca

    def cl_da(self):
        cl_da_duct = np.pi / 2 * self.zeta() * self.cl_a()
        print(f"cl_Da: {cl_da_duct}")
        return cl_da_duct

    def cl(self):
        k_prop = 0.2 * np.sqrt(np.abs(self.tc_prop))
        if self.pc == "off":
            norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
            norm_area = self.proj_area() / self.ref_area

            cl_duct = self.cl_da() * np.radians(self.inflow_angle())

            cl_duct_norm = cl_duct * norm_area * norm_speed
            return cl_duct, cl_duct_norm
        else:
            norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
            norm_area = self.proj_area() / self.ref_area

            cl_duct = (1 + k_prop) * self.cl_da() * np.radians(self.inflow_angle())

            cl_duct_norm = cl_duct * norm_area * norm_speed
            return cl_duct, cl_duct_norm

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), 't')
        fm = mach_correction(self.mach)
        ftc = 1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4
        coeff = airfoil_polar(f"support{self.duct_profile}.txt", float(0.0))
        cdmin = float(coeff[1] + coeff[2])
        # print(f"wet area: {self.wetted_area()}, ref: {self.ref_area}")
        cd0_duct = fm * ftc * cf * (self.wetted_area()/self.ref_area) # * (cdmin / 0.004) ** 4
        # cd0_duct = fm * ftc * cf * (cdmin / 0.004) ** 4
        return cd0_duct

    def cdi(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.proj_area() / self.ref_area

        oswald_duct = 2  # from vikesh

        cdi_duct = self.cl()[0] ** 2 / (oswald_duct * np.pi * self.aspect_ratio())

        cdi_norm = ((self.cl()[0] ** 2) / (oswald_duct * np.pi * self.aspect_ratio())) * norm_area * norm_speed
        return cdi_duct, cdi_norm

    def cd(self):
        cd_duct = self.cdi()[0] + self.cd0()
        return cd_duct

    def coefficient_norm(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.proj_area() / self.ref_area

        cd_norm = self.cd() * norm_speed * norm_area
        cl_norm = self.cl()[0] * norm_area * norm_speed

        return cl_norm, cd_norm

    """ ----------------------------------- determine output primes ---------------------------------------------- """
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


    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.proj_area() / self.ref_area

        cd_duct = self.cd() * norm_speed * norm_area
        return cd_duct

    def cl_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.proj_area() / self.ref_area

        cl_duct = self.cl()[0] * norm_speed * norm_area
        return cl_duct

    def cm(self):
        """ Aerodynamic moment of an annular wing based on Masqood"""
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.proj_area() / self.ref_area
        norm_chord = self.duct_chord / self.ref_chord

        alfa = np.radians(self.inflow_angle())
        ac = 0.3
        b = -0.786
        c = - 0.26

        p1 = 0.63
        p0 = -0.35

        kpcm = 6.25 * np.sin(self.aspect_ratio()/2)
        kvcm = np.pi / 3

        AR = 1.62

        xp = ac * self.aspect_ratio() ** b + c
        xe = p1 * self.aspect_ratio() + p0

        cm_duct = xp * kpcm * np.sin(alfa) * np.cos(alfa) + xe * kvcm * np.sin(alfa) ** 2
        cm_duct_norm = cm_duct * norm_chord * norm_speed * norm_area
        return cm_duct, cm_duct_norm

    def cn(self):
        """ return coefficient first, then return normalized coefficient"""
        n_area = self.proj_area() / self.ref_area
        n_velocity = self.inflow_velocity() ** 2 / self.v_inf ** 2
        alpha = np.radians(self.inflow_angle())

        cn_duct = self.cl() * np.cos(alpha) + self.cd() * np.sin(alpha)

        cn_norm = cn_duct * n_area * n_velocity
        return cn_duct, cn_norm

    def ct(self):
        """ return coefficient first, then return normalized coefficient"""
        n_area = self.proj_area() / self.ref_area
        n_velocity = self.inflow_velocity() ** 2 / self.v_inf ** 2
        alpha = np.radians(self.inflow_angle())

        ct_duct = self.cl() * np.sin(alpha) - self.cd() * np.cos(alpha)

        ct_norm = ct_duct * n_area * n_velocity
        return ct_duct, ct_norm

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

        a_c = -1 * np.cos(alpha)/(a_r * f4 +1)
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
        """ based on Torenbeek class II weight estimation"""
        kh = 1.05
        sh = self.duct_diameter * np.pi * self.duct_chord
        vd = ref.v_dive
        sweep = np.radians(0)

        m_duct = kh * sh * (62 * (sh ** 0.2 * vd) / (1000 * np.sqrt(np.cos(sweep))) - 2.5)
        return m_duct


""" Test section"""
"""
if __name__ == "__main__":

    a = np.linspace(0, 15, 31)
    cl = []
    a_ref = np.linspace(0, 15, 16)
    cl_ref = [0, 0.1, 0.203, 0.291, 0.357, 0.433, 0.511, 0.577, 0.702, 0.768, 0.844, 0.935, 1.067, 1.15, 1.2, 1.3]
    a_exp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
             29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
    cl_exp = [0.116, 0.198, 0.284, 0.369, 0.454, 0.54, 0.629, 0.719, 0.808, 0.898, 0.988, 1.083, 1.181, 1.28, 1.379,
              1.478, 1.577, 1.677, 1.776, 1.876, 1.976, 2.073, 2.17, 2.266, 2.363, 2.459, 2.446, 2.363, 2.281, 2.198,
              2.115, 2.093, 2.073, 2.053, 2.033, 2.012, 1.965, 1.917, 1.868, 1.82, 1.771, 1.725, 1.678, 1.632, 1.586,
              1.539, 1.492, 1.445, 1.398, 1.351, 1.304]
    cd_exp2 = [0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,1.02,1.04,1.06,1.08,1.1,1.12,1.14,1.16,1.18,1.2,1.22,1.24,1.26,1.28,1.3,1.32,1.34,1.36,1.38,1.4,1.42,1.44,1.46,1.48,1.5,1.52,1.54,1.56,1.58,1.6,1.62,1.64,1.66,1.68]
    cl_exp2 = [0.022,0.133,0.346,0.57,0.713,0.862,1.011,1.125,1.209,1.293,1.377,1.462,1.525,1.585,1.653,1.726,1.799,1.872,1.945,2.008,2.062,2.117,2.171,2.225,2.279,2.333,2.391,2.465,2.539,2.529,2.502,2.475,2.448,2.421,2.394,2.367,2.341,2.314,2.287,2.26,2.233,2.206,2.179,2.152,2.125,2.106,2.094,2.082,2.07,2.064,2.06,2.056,2.052,2.048,2.044,2.041,2.037,2.033,2.029,2.025,2.021,2.018,2.014,2,1.965,1.929,1.894,1.859,1.824,1.789,1.754,1.698,1.642,1.585,1.531,1.498,1.464,1.431,1.398,1.365,1.331]
    cl_the = []
    cd = []
    cd_0 = []
    cd_ref_appended = []
    cd_ref = [0.027, 0.028, 0.03, 0.033, 0.038, 0.041, 0.047, 0.054, 0.061, 0.071, 0.08, 0.087, 0.096, 0.109,
              0.125, 0.129]
    cd_ref= [0.025, 0.027, 0.03, 0.031, 0.037, 0.043, 0.046, 0.054, 0.064, 0.071, 0.084, 0.089, 0.095, 0.106, 0.125, 0.13]

    cm = []
    cm_ref = [0.002,-0.005,-0.009,-0.012,-0.015,-0.016,-0.016,-0.02,-0.022,-0.024,-0.024,-0.022,-0.02,-0.019,-0.017,-0.014]
    cm_dp = []
    ct_dp = []
    cn_dp = []

    cd_the = []
    for i in range(len(a)):
        wing = Duct(duct_diameter=config.duct_diameter,
                    duct_chord=config.duct_chord,
                    duct_profile=config.duct_airfoil,
                    alpha=a[i],
                    power_condition="on",
                    tc_prop=0.48,
                    v_inf=128,
                    mach=0.44,
                    ref_area=ref.s_w,
                    ref_chord=2.2345,
                    bem_input=[41420.85603250924, 26482.06279555917, -1.4475057750305072e-13, 0.8892292886261024, 0.3297344029147765, 0.3201225439053968, 5, 10],
                    va=723)
        al = np.radians(a[i])
        kp = 6.25 * np.sin(wing.aspect_ratio()/2)
        kv = np.pi / 3
        cl_theory = kp * np.sin(al) * np.cos(al)**2 + kv * np.cos(al) * np.sin(al) ** 2
        cl_the.append(cl_theory)
        cd_the.append(wing.cd0() + 0.06 * cl_theory ** 2)

        cl.append(wing.cl()[0])
        cd.append(wing.cd())
        cd_0.append(wing.cd0())
        cm.append(wing.cm()[0])
        cm_dp.append(wing.cm_dp())
        cn_dp.append(wing.cn_dp())
        ct_dp.append(wing.ct_dp())

    for j in range(len(a_ref)):
        cd_ref_appended.append(cd_ref[j] - 0.0231)  # correction for cd0

    alpha_avl, cl_avl, clff_avl, cd_avl, cdin_avl, cdff_avl = read_avl_output("AVL_RW_aeroproperties.txt")

    plt.figure('CL - alpha')
    plt.plot(a, cl, label=r'Model', color="tab:blue")
    plt.plot(a_ref, cl_ref, label=r'Experimental 1', color="tab:green", marker='o')
    plt.plot(a_exp, cl_exp, label=r'Experimental 2', color="tab:red", marker='o')
    plt.plot(a, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.plot(a, cn_dp, label="Nasa prediction", color="orange")
    plt.plot(alpha_avl, cl_avl, label="AVL", color="purple")
    plt.xlim([0, 15])
    plt.ylim([-0.1, 1.7])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Model', color="tab:blue")
    plt.plot(a_ref, cd_ref_appended, label=r'Experimental', color="tab:green", marker='o')
    plt.plot(a, cd_the, label='Theory', color="tab:green", linestyle="--")
    plt.plot(a, ct_dp, label="Nasa prediction", color="orange")
    plt.plot(alpha_avl, cd_avl, label="AVL", color="purple")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CM - alpha')
    plt.plot(a, cm, label=r'Model', color="tab:blue")
    plt.plot(a_ref, cm_ref, label=r'Experimental', color="tab:green", marker='o')  # for proper comparison change AR of the calculated values
    plt.plot(a, cm_dp, label="nasa", color="orange")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{M}$ [-]')
    plt.title(r'$C_{M}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD')
    plt.plot(cd, cl, label=r'Model', color="tab:blue")
    plt.plot(cd_ref_appended, cl_ref, label=r'Experimental 1', color="tab:green", marker='o')
    #plt.plot(cd_exp2, cl_exp2, label=r'Experimental 2', color="tab:red", marker='o')
    plt.plot(cd_the, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.plot(ct_dp, cn_dp, label="nasa", color="orange")
    plt.plot(cl_avl, cd_avl, label="AVL", color="purple")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.xlim([0, 0.25])
    plt.ylim([-0.1, 1.6])
    plt.title(r'$C_{L}$ vs. $C_{D}$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - CL^2')
    plt.plot([cl_val ** 2 for cl_val in cl], cd, label=r'Model', color="tab:blue")
    plt.plot([cl_val ** 2 for cl_val in cl_the], cd_the, label='Theory', color="tab:green", linestyle="--")
    plt.plot([cl_val ** 2 for cl_val in cl], cd_0, label=r'$C_{D0}$', color="tab:blue", linestyle="--")
    plt.plot([cl_avl ** 2 for cl_avl in cl_avl], cd_avl, label="AVL", color="purple")
    plt.xlabel(r'$C_{L}^{2}$ [-]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $C_{L}^{2}$ - Duct')
    plt.xlim([0, 0.5])
    plt.ylim([0, 0.04])
    plt.legend()
    plt.grid(True)
    plt.show()

    print(f"inflow vel: {wing.inflow_velocity()}")
    print(f"inflow ang: {wing.inflow_angle()}")
    print(f"area: {wing.proj_area()}, wetted area: {wing.wetted_area()}")
    print(f"aspect ratio: {wing.aspect_ratio()}, t_c: {wing.t_c()}")
    print(f"cd0: {wing.cd0()}, cdi: {wing.cdi()}, cdprime: {wing.cd_prime()}")
    print(f"cl: {wing.cl()}, cl_prime: {wing.cl_prime()}")
    print(f"weight: {wing.weight()}") """
