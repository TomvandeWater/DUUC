import numpy as np
from analysis_modules.factors import skin_friction, mach_correction
from analysis_modules.aerodynamic import reynolds
from analysis_modules.ISA import air_density_isa
import config
import data.atr_reference as ref
import matplotlib.pyplot as plt


class Nacelle:

    def __init__(self, geometry, conditions, reference, propulsor_type: str, power_condition: str, v_after_prop: float,
                 a_after_prop: float):
        super().__init__()

        self.nacelle_length = geometry[0]
        self.nacelle_diameter = geometry[1]

        self.propulsor_type = propulsor_type
        self.pc = power_condition

        self.v_inf = conditions[0]
        self.alpha = conditions[1]
        self.altitude = conditions[2]
        self.mach = conditions[3]
        self.density = air_density_isa(self.altitude)[0]

        self.ref_area = reference[0]
        self.ref_chord = reference[1]

        self.a_after_prop = a_after_prop
        self.v_after_prop = v_after_prop

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

        angle_rad = np.radians(inflow_nacelle)
        return inflow_nacelle, angle_rad

    def reynolds_number(self):
        re_nac = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.nacelle_length)
        return re_nac

    """ -------------------------------------- geometric properties ---------------------------------------------- """
    def area_wetted(self):
        """ only the side area is assumed no closing sides"""
        area_cylinder = np.pi * self.nacelle_diameter * self.nacelle_length

        """ half a sphere is used to close of the rear end"""
        area_rear = 0.5 * np.pi * self.nacelle_diameter ** 2

        area_nacelle = area_rear + area_cylinder
        return area_nacelle

    def area(self):
        """ projected area on the bottom """
        area_nacelle = self.nacelle_length * self.nacelle_diameter
        return area_nacelle

    """ ------------------------------ Ratios used for normalization --------------------------------------------- """
    def area_ratio(self):
        ar_nacelle = self.area() / self.ref_area
        return ar_nacelle

    def area_ratio_wet(self):
        ar_w_nacelle = self.area_wetted() / self.ref_area
        return ar_w_nacelle

    def velocity_ratio(self):
        v_ratio = self.inflow_velocity() ** 2 / self.v_inf ** 2
        return v_ratio

    """ ---------------------------------------- determine coefficients ------------------------------------------ """
    def cd0(self):

        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)
        l_d = self.nacelle_length / self.nacelle_diameter
        f_nac = 1 + 60 / l_d ** 3 + 0.0025 * l_d

        cd0_nacelle = cf * fm * f_nac
        cd0_nacelle_ar = cf * fm * f_nac * self.area_ratio_wet()
        return cd0_nacelle, cd0_nacelle_ar

    def cd(self):
        """ drag model based on flow over slender inclined bodies of revolution """
        alpha = self.inflow_angle()[1]
        area = self.ref_area
        area_p = self.area()
        surface_b = self.area()
        cd_90 = 0.19

        cd_nac = (self.cd0()[1] * np.cos(alpha) ** 3 + surface_b / area * np.sin(2 * alpha) * np.sin(alpha / 2)
                  + cd_90 * area_p / area * np.sin(alpha) ** 3)
        cd_norm = cd_nac * self.velocity_ratio()

        # create vector for plots in GUI
        a_vector = np.linspace(-5, 15, 21)
        cd_vector = []
        for i in range(len(a_vector)):
            a_rad = np.radians(a_vector[i])
            cd_vector.append(self.cd0()[1] * np.cos(a_rad) ** 3 + surface_b / area * np.sin(2 * a_rad)
                             * np.sin(a_rad / 2) + cd_90 * area_p / area * np.sin(a_rad) ** 3)

        return cd_nac, cd_norm, cd_vector

    def cl(self):
        """ lift model based on flow over slender inclined bodies of revolution """
        alpha = self.inflow_angle()[1]
        area = self.ref_area
        area_p = self.area()
        surface_b = self.area()
        cd_90 = 0.19

        cl_nac = (surface_b / area * np.sin(2 * alpha) * np.cos(alpha / 2) + cd_90 * area_p / area * np.sin(alpha) ** 2
                  * np.cos(alpha))
        cl_norm = cl_nac * self.area_ratio_wet()

        # vectors for plots in GUI
        a_vector = np.linspace(-5, 15, 21)
        cl_vector = []
        for i in range(len(a_vector)):
            a_rad = np.radians(a_vector[i])
            cl_vector.append(surface_b / area * np.sin(2 * a_rad) * np.cos(a_rad / 2) + cd_90 * area_p / area
                             * np.sin(a_rad) ** 2 * np.cos(a_rad))

        return cl_nac, cl_norm, cl_vector

    def cm(self):
        """ pitching moment based on flow over slender inclined bodies of revolution """
        alpha = self.inflow_angle()[1]
        area = self.ref_area
        area_p = self.area()
        surface_b = self.area()
        chord = self.ref_chord
        cd_90 = 0.19
        length = self.nacelle_length
        xm = 0
        q = np.pi * (self.nacelle_diameter / 2) ** 2 * self.nacelle_length
        xa90 = self.nacelle_length / 2

        cm_nac = ((q - surface_b * (length - xm) / (area * chord)) * np.sin(2 * alpha) * np.cos(alpha / 2) + cd_90
                  * (area_p / area) * ((xm - xa90) / chord) * np.sin(alpha) ** 2)
        cm_norm = cm_nac * self.velocity_ratio()

        cm_vector = []
        a_vect = np.linspace(-5, 15, 21)
        for i in range(len(a_vect)):
            cm_vector.append(((q - surface_b * (length - xm) / (area * chord)) * np.sin(2 * alpha) * np.cos(alpha / 2)
                              + cd_90 * (area_p / area) * ((xm - xa90) / chord) * np.sin(alpha) ** 2))

        return cm_nac, cm_norm, cm_vector

    def cn(self):
        alpha = self.inflow_angle()[1]

        cn_nacelle = self.cl()[0] * np.cos(alpha) + self.cd()[0] * np.sin(alpha)

        cn_norm = cn_nacelle * self.area_ratio() * self.velocity_ratio()
        return cn_nacelle, cn_norm

    def ct(self):
        alpha = self.inflow_angle()[1]

        ct_nacelle = self.cl()[0] * np.sin(alpha) - self.cd()[0] * np.cos(alpha)

        ct_norm = ct_nacelle * self.area_ratio() * self.velocity_ratio()

        return ct_nacelle, ct_norm

    """ ------------------------------------------- FORCES --------------------------------------------------------- """
    def lift_force(self):
        lift_nacelle = self.cl()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return lift_nacelle

    def drag_force(self):
        drag_nacelle = self.cd()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return drag_nacelle

    def moment_force(self):
        moment_nacelle = self.cm()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return moment_nacelle

    """ -------------------------------------- The weight is depending on the propulsion type ---------------------- """
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
