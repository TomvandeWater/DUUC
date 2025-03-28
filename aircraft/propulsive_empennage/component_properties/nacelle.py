import numpy as np
from analysis_modules.factors import skin_friction, mach_correction
from analysis_modules.aerodynamic import reynolds
from analysis_modules.ISA import air_density_isa
import config
import data.atr_reference as ref
import matplotlib.pyplot as plt


class Nacelle:

    def __init__(self, nacelle_length: float, nacelle_diameter: float,
                 propulsor_type: str, power_condition: str, v_after_prop: float,
                 alpha: float, ref_area: float, v_inf: float, mach: float, a_after_prop: float, altitude: float):
        super().__init__()
        self.nacelle_length = nacelle_length
        self.nacelle_diameter = nacelle_diameter
        self.propulsor_type = propulsor_type
        self.pc = power_condition
        self.v_after_prop = v_after_prop
        self.alpha = alpha
        self.ref_area = ref_area
        self.v_inf = v_inf
        self.mach = mach
        self.a_after_prop = a_after_prop
        self.altitude = altitude
        self.density = air_density_isa(self.altitude)[0]

    """ -------------------------------- inflow properties ------------------------------------------------------ """
    def inflow_velocity(self):
        if self.pc == "off":
            u_nacelle = self.v_inf
            return u_nacelle
        else:
            """ Effective velocity after the propeller blade influenced by the support strut """
            v1 = np.cos(np.radians(self.a_after_prop)) * self.v_after_prop
            v2 = np.cos(np.radians(self.a_after_prop * 0.5)) * self.v_after_prop

            u_nacelle = (v1 + v2) / 2
            return u_nacelle

    def inflow_angle(self):
        """ This is 0.75 * the angle of the flow after the propeller"""
        inflow_nacelle = self.a_after_prop * 0.75

        # inflow_nacelle = self.alpha  # uncomment this for alpha variation curves.
        inflow_nacelle = self.alpha
        return inflow_nacelle

    def reynolds_number(self):
        re_nac = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.nacelle_length)
        return re_nac

    """ -------------------------------------- geometric properties ---------------------------------------------- """
    def wet_area(self):
        """ only the side area is assumed no closing sides"""
        area_cylinder = np.pi * self.nacelle_diameter * self.nacelle_length

        """ half a sphere is used to close of the rear end"""
        area_rear = 0.5 * np.pi * self.nacelle_diameter ** 2

        area_nacelle = area_rear + area_cylinder
        return area_nacelle

    def area(self):
        """ projected area"""
        area_nacelle = self.nacelle_length * self.nacelle_diameter
        return area_nacelle

    """ ---------------------------------------- determine coefficients ------------------------------------------ """
    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)
        l_d = self.nacelle_length / self.nacelle_diameter
        f_nac = 1 + 60 / l_d ** 3 + 0.0025 * l_d
        area_ratio = (self.nacelle_length * self.nacelle_diameter) / self.ref_area

        cd0_nacelle = cf * fm * f_nac * area_ratio

        return cd0_nacelle

    @staticmethod
    def cdi():
        """ assume nacelle does not produce lift"""
        cdi_nac = 0
        return cdi_nac

    def cd(self):
        cd_nac = self.cd0() + self.cdi()

        alpha = np.radians(self.inflow_angle())
        A = self.ref_area
        Ap = self.area()
        Sb = self.area()
        cd_90 = 0.19

        cd_2 = (self.cd0() * np.cos(alpha) ** 3 + Sb / A * np.sin(2 * alpha) * np.sin(alpha / 2)
                + cd_90 * Ap / A * np.sin(alpha) ** 3)
        return cd_nac, cd_2

    def cl(self):
        cl_nac = 0

        alpha = np.radians(self.inflow_angle())
        A = self.ref_area
        Ap = self.area()
        Sb = self.area()
        cd_90 = 0.19

        cl_2 = Sb / A * np.sin(2 * alpha) * np.cos(alpha / 2) + cd_90 * Ap / A * np.sin(alpha) ** 2 * np.cos(alpha)
        return cl_nac, cl_2

    def coefficients_norm(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.wet_area() / self.ref_area

        cd_norm = self.cd()[1] * norm_speed * norm_area
        cl_norm = 0

        return cl_norm, cd_norm

    """ ------------------------------------------- prime outputs ------------------------------------------------ """
    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.wet_area() / self.ref_area

        cd_nacelle = self.cd()[1] * norm_speed
        return cd_nacelle

    def cl_prime(self):
        """ assume no lift is produced by the nacelle """
        cl_nacelle = self.cl()[1]
        return cl_nacelle

    def cm(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.wet_area() / self.ref_area
        """ based on Htail mounted propeller research by van Arnhem 2019 (Aerodynamic Performance of an Aircraft
        Equipped with Horizontal Tail Mounted Propellers)"""
        a = 0.033
        b = -0.013

        cm_nac = a * self.inflow_angle() + b
        cm_norm = cm_nac * norm_area * norm_speed
        return cm_nac, cm_norm

    def cn(self):
        n_area = self.area() / self.ref_area
        n_velocity = self.inflow_velocity() ** 2 / self.v_inf ** 2

        alpha = np.radians(self.inflow_angle())

        cn_nacelle = self.cl()[1] * np.cos(alpha) + self.cd()[1] * np.sin(alpha)

        cn_norm = cn_nacelle * n_area * n_velocity
        return cn_nacelle, cn_norm

    def ct(self):
        n_area = self.area() / self.ref_area
        n_velocity = self.inflow_velocity() ** 2 / self.v_inf ** 2

        alpha = np.radians(self.inflow_angle())

        ct_nacelle = self.cl()[1] * np.sin(alpha) - self.cd()[1] * np.cos(alpha)

        ct_norm = ct_nacelle * n_area * n_velocity

        return ct_nacelle, ct_norm

    """ -------------------------------------- The weight is depending on the propulsor type --------------------- """

    def weight(self):
        if self.propulsor_type == 'conventional':
            """ based on Torenbeek Class II weight estimation"""
            # calculated for take off conditions
            p = 4102000
            eta = 0.64545
            v = 54.12

            m_nacelle = 0.0458 * (p * eta) / (v * 9.81)

            return m_nacelle

        else:
            """ based on Torenbeek Class II weight estimation"""
            # calculated for take off conditions
            p = 4102000
            eta = 0.64545
            v = 54.12

            m_nacelle = 0.0458 * (p * eta) / (v * 9.81)

            return m_nacelle


""" Test section """
"""
if __name__ == "__main__":
    a = np.linspace(0, 15, 31)
    cl = []
    cd = []
    cm = []
    cd_2 = []
    cl_2 = []

    for i in range(len(a)):
        nacelle = Nacelle(nacelle_length=config.nacelle_length,
                          nacelle_diameter=config.nacelle_diameter,
                          propulsor_type="conventional",
                          power_condition="on",
                          alpha=a[i],
                          ref_area=ref.s_w,
                          v_inf=128,
                          mach=0.44,
                          a_after_prop=20,
                          v_after_prop=125)

        cl.append(nacelle.cl()[0])
        cd.append(nacelle.cd()[0])
        cm.append(nacelle.cm()[1])
        cl_2.append(nacelle.cl()[1])
        cd_2.append(nacelle.cd()[1])

    plt.figure('CL - alpha')
    plt.plot(a, cl, label=r'Prediction model', color="tab:blue")
    plt.plot(a, cl_2, label=r'ILR Reference', color="tab:blue", linestyle="dashed")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Nacelle')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Prediction model', color="tab:blue")
    plt.plot(a, cd_2, label=r'ILR Reference', color="tab:blue", linestyle="dashed")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'Nacelle')
    plt.legend()
    plt.grid(True)

    plt.figure('CM - alpha')
    plt.plot(a, cm, label=r'Prediction model', color="tab:blue")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{M}$ [-]')
    plt.title(r'Nacelle')
    plt.legend()
    plt.grid(True)
    plt.show()

    print(f"inflow vel: {nacelle.inflow_velocity()}")
    print(f"inflow ang: {nacelle.inflow_angle()}")
    print(f"cd0: {nacelle.cd0()}")
    print(f"cd prime: {nacelle.cd_prime():.3f}")
"""
