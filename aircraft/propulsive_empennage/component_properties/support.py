import numpy as np
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa
from data.read_data import airfoil_polar
import config
from analysis_modules.factors import skin_friction, mach_correction
import matplotlib.pyplot as plt
import data.atr_reference as ref


class SupportStrut:
    def __init__(self, support_length: float, support_chord: float, support_profile: str,
                 cant_angle: float, power_condition: str, v_after_prop: float,
                 alpha: float, tc_prop: float, cn_prop: float, ref_area: float, v_inf: float,
                 va_inlet: float, prop_diameter: float, a_after_prop: float, m_supported: float,
                 ref_chord: float, altitude: float):
        super().__init__()
        self.support_length = support_length
        self.support_chord = support_chord
        self.support_profile = support_profile
        self.cant_angle = cant_angle
        self.pc = power_condition
        self.v_after_prop = v_after_prop
        self.alpha = alpha
        self.tc_prop = tc_prop
        self.cn_prop = cn_prop
        self.ref_area = ref_area
        self.v_inf = v_inf
        self.va_inlet = va_inlet
        self.prop_diameter = prop_diameter
        self.a_after_prop = a_after_prop
        self.m_supported = m_supported
        self.ref_chord = ref_chord
        self.altitude = altitude
        self.density = air_density_isa(self.altitude)[0]

    """ ------------------ inflow properties ----------------------------------------------------------  """
    def inflow_velocity(self):
        if self.pc == "off":
            u_support = self.va_inlet / (np.pi / 4 * self.prop_diameter ** 2)
            return u_support
        else:
            """ Effective velocity after the propeller"""
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
            cant = np.radians(self.cant_angle)

            inflow_support = self.a_after_prop * np.cos(cant)
            return inflow_support

    def reynolds_number(self):
        re_sup = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.support_chord)
        return re_sup

    """  ------------------------------------- geometric properties -------------------------------------------- """
    def area(self):
        """The area is based on a rectangle, note that the support goes through the nacelle in this
        model. The area is hence overestimated."""
        s_support = self.support_chord * self.support_length
        return s_support

    def area_wetted(self):
        s_wet_support = 2 * (1 + 0.5 * self.t_c()) * self.support_length * self.support_chord
        return s_wet_support

    def t_c(self):
        num_list = [int(digit) for digit in self.support_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def aspect_ratio(self):
        ar_support = self.support_length ** 2 / self.area()
        return ar_support

    """ -------------------------------------- coefficients for support forces ------------------------------------ """
    def cl_da(self):
        cl0 = airfoil_polar(f"support{self.support_profile}.txt", float(0.0))
        cl0_val = float(cl0[0])
        cl5 = airfoil_polar(f"support{self.support_profile}.txt", float(10.0))
        cl5_val = float(cl5[0])
        cl_da_support = (cl5_val - cl0_val) / 10
        return cl_da_support

    def cl(self):
        cl_support = self.cl_da() * self.inflow_angle()
        return cl_support

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(0.44)

        ftc = (1 + 1.2 * self.t_c() + 60 * self.t_c() ** 4)

        cd0_support = cf * fm * ftc * (self.area_wetted() / self.ref_area)

        return cd0_support

    def cdi(self):
        cdi_support = self.cl()[0] ** 2 / (0.95 * np.pi * self.aspect_ratio())
        return cdi_support

    def cd(self):
        cd_polar = airfoil_polar(f"support{self.support_profile}.txt", float(self.inflow_angle()))
        cd_val = float(cd_polar[1]) + self.cd0()

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

    def coefficient_norm(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_norm = self.cd() * norm_speed * norm_area
        cl_norm = self.cl() * norm_speed * norm_area
        return cl_norm, cd_norm

    """ ------------------------------- determine output primes ------------------------------------------------- """
    def cl_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cl_support = self.cl() * norm_area * norm_speed
        return cl_support

    def cd_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_support = self.cd() * norm_area * norm_speed
        return cd_support

    def cm(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_chord = self.support_chord / self.ref_chord

        cm_strut = airfoil_polar(f"support{self.support_profile}.txt", float(self.inflow_angle()))
        cm_strut = float(cm_strut[2])
        cm_norm = cm_strut * norm_chord * norm_speed * norm_area
        return cm_strut, cm_norm

    def cn(self):
        n_area = self.area() / self.ref_area
        n_velocity = self.inflow_velocity() ** 2 / self.v_inf ** 2

        alpha = np.radians(self.inflow_angle())

        cn_support = self.cl() * np.cos(alpha) + self.cd() * np.sin(alpha)

        cn_norm = cn_support * n_area * n_velocity
        return cn_support, cn_norm

    def ct(self):
        n_area = self.area() / self.ref_area
        n_velocity = self.inflow_velocity() ** 2 / self.v_inf ** 2

        alpha = np.radians(self.inflow_angle())

        ct_support = self.cl() * np.sin(alpha) - self.cd() * np.cos(alpha)

        ct_norm = ct_support * n_velocity * n_area

        return ct_support, ct_norm

    """ ----------------------------------- Weight estimation ------------------------------------------------- """
    @staticmethod
    def weight():
        """ Weight is calculated in the empennage_assembly_PE.py file"""
        m_support = None
        return m_support


""" Test section """

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
                               ref_area=ref.s_w,
                               v_inf=128,
                               a_after_prop=a[i],
                               m_supported=1000,
                               prop_diameter=3.6,
                               va_inlet=250,
                               ref_chord=2.2345)
        """
        polar = airfoil_polar(f"support0012.txt", float(a[i]))
        cd_val = float(polar[1])
        cl_val = float(polar[0])

        al = np.radians(a[i])
        cl_theory = np.pi * 2 * al
        cl_the.append(cl_theory)
        cd_the.append(support.cd0())

        cl.append(support.cl())
        cd.append(support.cd())
        cd_ref.append(cd_val[0])
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
    """
    print(f"inflow vel: {support.inflow_velocity()}")
    print(f"inflow ang: {support.inflow_angle()}")
    print(f"cd0: {support.cd0()}")
    print(f"cd interference: {support.cd_interference():.5f}")
    print(f"cd prime: {support.cd_prime():.5f}")
    print(f"cl: {support.cl():.5f}")
    print(f"cl prime: {support.cl_prime():.5f}")
    print(f"wet area: {support.area_wetted()}")
