import numpy as np
from analysis_modules.factors import mach_correction
from analysis_modules.factors import skin_friction
from data.read_data import airfoil_polar
import data.atr_reference as ref
import config
import matplotlib.pyplot as plt


class Duct:
    def __init__(self, duct_diameter: float, duct_chord: float, duct_profile: str,
                 alpha: float, re_duct: float, power_condition: str, u_mom: float,
                 tc_prop: float, v_inf: float, mach: float, ref_area: float):
        super().__init__()
        self.duct_diameter = duct_diameter
        self.duct_chord = duct_chord
        self.duct_profile = duct_profile
        self.alpha = alpha
        self.re_duct = re_duct
        self.pc = power_condition
        self.u_mom = u_mom
        self.tc_prop = tc_prop
        self.v_inf = v_inf
        self.mach = mach
        self.ref_area = ref_area

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
        return cl_da_duct

    def cl(self):
        k_prop = 0.2 * np.sqrt(self.tc_prop)
        if self.pc == "off":
            cl_duct = self.cl_da() * np.radians(self.inflow_angle())
            return cl_duct
        else:
            cl_duct = (1 + k_prop) * self.cl_da() * np.radians(self.inflow_angle())
            return cl_duct

    def cd0(self):
        cf = skin_friction(self.re_duct, 'T')
        fm = mach_correction(self.mach)
        ftc = 1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4
        coeff = airfoil_polar(f"support{self.duct_profile}.txt", float(0.0))
        cdmin = float(coeff[1] + coeff[2])

        cd0_duct = fm * ftc * cf * (self.wetted_area()/self.proj_area()) * (cdmin / 0.004) ** 4
        return cd0_duct

    def cdi(self):
        cdi_duct = self.cl() ** 2 / (2 * np.pi * self.aspect_ratio())
        return cdi_duct

    def cd(self):
        cd_duct = self.cdi() + self.cd0()
        return cd_duct

    """ ----------------------------------- determine output primes ---------------------------------------------- """
    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.proj_area() / self.ref_area

        cd_duct = self.cd() * norm_speed * norm_area
        return cd_duct

    def cl_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.proj_area() / self.ref_area

        cl_duct = self.cl() * norm_speed * norm_area
        return cl_duct

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
    cd_ref = [0.027, 0.028, 0.03, 0.033, 0.038, 0.041, 0.047, 0.054, 0.061, 0.071, 0.08, 0.087, 0.096, 0.109,
              0.125, 0.129]
    cd_the = []
    for i in range(len(a)):
        wing = Duct(duct_diameter=config.duct_diameter,
                    duct_chord=config.duct_chord,
                    duct_profile=config.duct_airfoil,
                    alpha=a[i],
                    re_duct=8422274,
                    power_condition="on",
                    u_mom=131,
                    tc_prop=0.48,
                    v_inf=128,
                    mach=0.44,
                    ref_area=ref.s_w)
        al = np.radians(a[i])
        kp = 6.25 * np.sin(wing.aspect_ratio()/2)
        kv = np.pi / 3
        cl_theory = kp * np.sin(al) * np.cos(al)**2 + kv * np.cos(al) * np.sin(al) ** 2
        cl_the.append(cl_theory)
        cd_the.append(wing.cd0() + 0.06 * cl_theory ** 2)

        cl.append(wing.cl())
        cd.append(wing.cd())

    plt.figure('CL - alpha')
    plt.plot(a, cl, label=r'Model', color="tab:blue")
    plt.plot(a_ref, cl_ref, label=r'Experimental 1', color="tab:green", marker='o')
    plt.plot(a_exp, cl_exp, label=r'Experimental 2', color="tab:red", marker='o')
    plt.plot(a, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlim([0, 15])
    plt.ylim([-0.1, 1.7])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Model', color="tab:blue")
    plt.plot(a_ref, cd_ref, label=r'Experimental', color="tab:green", marker='o')
    plt.plot(a, cd_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD')
    plt.plot(cd, cl, label=r'Model', color="tab:blue")
    plt.plot(cd_ref, cl_ref, label=r'Experimental 1', color="tab:green", marker='o')
    plt.plot(cd_exp2, cl_exp2, label=r'Experimental 2', color="tab:red", marker='o')
    plt.plot(cd_the, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.xlim([0, 0.25])
    plt.ylim([-0.1, 1.6])
    plt.title(r'$C_{L}$ vs. $C_{D}$ - Duct')
    plt.legend()
    plt.grid(True)
    plt.show()

    print(f"inflow vel: {wing.inflow_velocity()}")
    print(f"inflow ang: {wing.inflow_angle()}")
    print(f"area: {wing.proj_area()}, wetted area: {wing.wetted_area()}")
    print(f"aspect ratio: {wing.aspect_ratio()}, t_c: {wing.t_c()}")
    print(f"cd0: {wing.cd0()}, cdi: {wing.cdi()}, cdprime: {wing.cd_prime()}")
    print(f"cl: {wing.cl()}, cl_prime: {wing.cl_prime()}")
    print(f"weight: {wing.weight()}")
