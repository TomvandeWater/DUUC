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
    plt.plot(a_ref, cl_ref, label=r'Experimental', color="tab:green", marker='o')
    plt.plot(a, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'CL - $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Model', color="tab:blue")
    plt.plot(a_ref, cd_ref, label=r'Experimental', color="tab:green", marker='o')
    plt.plot(a, cd_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'CD - $\alpha$ - Duct')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD')
    plt.plot(cd, cl, label=r'Model', color="tab:blue")
    plt.plot(cd_ref, cl_ref, label=r'Experimental', color="tab:green", marker='o')
    plt.plot(cd_the, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ - $C_{D}$')
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
