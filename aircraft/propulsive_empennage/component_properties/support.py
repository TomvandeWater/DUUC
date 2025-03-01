import numpy as np
from analysis_modules.aerodynamic import drag_interference
from data.read_data import airfoil_polar
import config
from analysis_modules.factors import skin_friction, mach_correction
import matplotlib.pyplot as plt


class SupportStrut:
    def __init__(self, support_length: float, support_chord: float, support_profile: str,
                 cant_angle: float, power_condition: str, v_after_prop: float, u_mom: float,
                 alpha: float, tc_prop: float, cn_prop: float, ref_area: float, v_inf: float,
                 va_inlet: float, prop_diameter: float, a_after_prop: float, m_supported: float):
        super().__init__()
        self.support_length = support_length
        self.support_chord = support_chord
        self.support_profile = support_profile
        self.cant_angle = cant_angle
        self.pc = power_condition
        self.v_after_prop = v_after_prop
        self.u_mom = u_mom
        self.alpha = alpha
        self.tc_prop = tc_prop
        self.cn_prop = cn_prop
        self.ref_area = ref_area
        self.v_inf = v_inf
        self.va_inlet = va_inlet
        self.prop_diameter = prop_diameter
        self.a_after_prop = a_after_prop
        self.m_supported = m_supported

    """ ------------------ inflow properties ----------------------------------------------------------  """
    def inflow_velocity(self):
        if self.pc == "off":
            u_support = self.va_inlet / (np.pi / 4 * self.prop_diameter ** 2)
            return u_support
        else:
            u_support = self.v_after_prop
            return u_support

    def inflow_angle(self):
        thrust_c = self.tc_prop
        c_norm = self.cn_prop
        de_da_prop = ((thrust_c / (4 + 8/7 * thrust_c)) + (c_norm * (1 + 1.3 * thrust_c) ** 0.5)
                      / (4 + 2 * thrust_c))

        if self.pc == "off":
            inflow_support = 0
            return inflow_support
        else:
            inflow_support = self.a_after_prop
            return inflow_support

    """  ------------------------------------- geometric properties -------------------------------------------- """
    def area(self):
        """The area is based on a rectangle, note that the support goes through the nacelle in this
        model. The area is hence overestimated."""
        s_support = self.support_chord * self.support_length
        return s_support

    def t_c(self):
        num_list = [int(digit) for digit in self.support_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    """ -------------------------------------- coefficients for support forces ------------------------------------ """
    def cl_da(self):
        cl0 = airfoil_polar(f"support{self.support_profile}.txt", float(0.0))
        cl0_val = float(cl0[0])
        cl5 = airfoil_polar(f"support{self.support_profile}.txt", float(10.0))
        cl5_val = float(cl5[0])
        cl_da_support = (cl5_val - cl0_val) / 10
        # print(f"dcl da: {cl_da_support}")
        return cl_da_support

    def cl(self):
        cl_support = self.cl_da() * self.inflow_angle()
        # print(f"inflow angle support: {self.inflow_angle()}")
        # print(f"support cl: {cl_support}")
        return cl_support

    def cd0(self):
        cf = skin_friction(8422274, "t")
        fm = mach_correction(0.44)

        ftc = (1 + 1.2 * self.t_c() + 60 * self.t_c() ** 4)

        cd0_support = cf * fm * ftc

        return cd0_support

    def cd(self):
        cd = airfoil_polar(f"support{self.support_profile}.txt", float(self.inflow_angle()))
        cd_val = float(cd[1]) + self.cd0()
        return cd_val

    def cd_interference(self):
        norm_area = (self.t_c() * self.support_chord ** 2) / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_support_nacelle = (2 * drag_interference(self.t_c(), "plane")
                              * norm_area * norm_speed)  # multiplied by 2 for 2 nac-supp inter

        cd_support_duct = (2 * drag_interference(self.t_c(), "t-junction")
                           * norm_area * norm_speed)  # multiplied by 2 for 2 supp-duct inter

        cd_int_support = cd_support_nacelle + cd_support_duct
        return cd_int_support

    """ ------------------------------- determine output primes ------------------------------------------------- """
    def cl_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cant = np.radians(self.cant_angle)

        cl_cl = (self.cl() * np.cos(cant) * np.cos(np.radians(self.inflow_angle())) * norm_speed
                 * norm_area)

        cl_cd = (self.cd() * np.sin(np.radians(self.inflow_angle())) * norm_speed
                 * norm_area)

        cl_support = cl_cl + cl_cd
        return cl_support

    def cd_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cant = np.radians(self.cant_angle)

        cd_cd = (self.cd() * np.cos(np.radians(self.inflow_angle())) * norm_speed
                 * norm_area)

        cd_cl = (self.cl() * np.cos(cant) * np.sin(np.radians(self.inflow_angle())) * norm_speed
                 * norm_area)

        cd_support = cd_cd + cd_cl
        return cd_support

    def cm(self):
        cm_strut = - self.cl() * 0.25
        return cm_strut

    """ ----------------------------------- Weight estimation ------------------------------------------------- """
    def weight(self):
        m_support = 0.1 * self.m_supported
        return m_support


""" Test section """
"""
if __name__ == "__main__":

    a = np.linspace(1, 2, 3)
    cl = []
    #a_ref = np.linspace(0, 15, 16)
    cl_ref = []
    cl_the = []
    cd = []
    cd_ref = []
    cd_the = []
    for i in range(len(a)):
        support = SupportStrut(support_length=config.support_length,
                               support_chord=config.support_chord,
                               support_profile=config.support_airfoil,
                               cant_angle=config.cant_angle,
                               power_condition="on",
                               v_after_prop=130,
                               u_mom=110,
                               alpha=0,
                               tc_prop=0.37,
                               cn_prop=0.09,
                               ref_area=config.duct_chord * config.duct_diameter,
                               v_inf=128,
                               a_after_prop=a[i],
                               m_supported=1000,
                               prop_diameter=3.6,
                               va_inlet=250)
        polar = airfoil_polar(f"support0012.txt", float(a[i]))
        cd_val = float(polar[1] + polar[2])
        cl_val = float(polar[0])

        al = np.radians(a[i])
        cl_theory = np.pi * 2 * al
        cl_the.append(cl_theory)
        cd_the.append(support.cd0())

        cl.append(support.cl())
        cd.append(support.cd())
        cd_ref.append(cd_val)
        cl_ref.append(cl_val)

    plt.figure('CL - alpha')
    plt.plot(a, cl, label=r'Model', color="tab:blue")
    plt.plot(a, cl_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(a, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Support')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Model', color="tab:blue")
    plt.plot(a, cd_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(a, cd_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'Support')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD')
    plt.plot(cd, cl, label=r'Model', color="tab:blue")
    plt.plot(cd_ref, cl_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(cd_the, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Support')
    plt.legend()
    plt.grid(True)
    plt.show()

    print(f"inflow vel: {support.inflow_velocity()}")
    print(f"inflow ang: {support.inflow_angle()}")
    print(f"cd: {support.cd():.5f}")
    print(f"cd interference: {support.cd_interference():.5f}")
    print(f"cd prime: {support.cd_prime():.5f}")
    print(f"cl: {support.cl():.5f}")
    print(f"cl prime: {support.cl_prime():.5f}") """
