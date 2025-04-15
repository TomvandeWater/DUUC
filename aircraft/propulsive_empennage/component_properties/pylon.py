import numpy as np
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa
import matplotlib.pyplot as plt
from data.read_data import airfoil_polar
import data.atr_reference as ref
import config
from analysis_modules.factors import skin_friction, mach_correction


class Pylon:
    """ The reference point for connecting the pylon to the fuselage is on y=0 and in at half of the chord length of
    the airfoil. The pylon is completely outside the duct. """
    def __init__(self, geometry, conditions, reference, power_condition: str):
        super().__init__()
        self.pylon_length = geometry[0]
        self.pylon_chord = geometry[1]
        self.pylon_airfoil = geometry[2]
        self.cant_angle = geometry[3]

        self.v_inf = conditions[0]
        self.alpha = conditions[1]
        self.altitude = conditions[2]
        self.mach = conditions[3]
        self.density = air_density_isa(self.altitude)[0]
        self.beta = conditions[6]

        self.ref_area = reference[0]
        self.ref_chord = reference[1]

        self.pc = power_condition

    """ The inflow speed for the pylon is affected by the outside of the duct and also affected by the downwash of the 
    wing. The downwash of the wing and the duct affect the inlfow angle of on the pylon"""
    def inflow_velocity(self):
        """ Inflow velocity of the pylon """
        e = (2 * ref.cl_wing) / (np.pi * ref.ar_w)

        if self.pc == "off":
            u_pylon = self.v_inf * (1 - (e ** 2) / 2)
            return u_pylon
        else:
            u_pylon = self.v_inf * (1 - (e ** 2) / 2)
            return u_pylon

    def inflow_angle(self):
        """ inflow angle already determined in general propulsive empennage class is only corrected here for the
        cant angle of the pylon -> asuming sweep = 0 """
        cant = np.radians(self.cant_angle)

        inflow_angle = self.alpha + - cant * np.cos(0)
        angle_rad = np.radians(inflow_angle)
        return inflow_angle, angle_rad

    def reynolds_number(self):
        """ local reynolds number of the pylon"""
        re_pylon = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.pylon_chord)
        return re_pylon

    """ ------------------------- geometric properties ----------------------------------------------------------- """
    def area(self):
        """ flat surface area of the pylon"""
        area_pylon = self.pylon_length * self.pylon_chord
        return area_pylon

    def t_c(self):
        """ assume NACA 4-series airfoil """
        num_list = [int(digit) for digit in self.pylon_airfoil]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def area_wetted(self):
        """ untapered, unswept wing definition by Torenbeek """
        s_wet_pylon = 2 * (1 + 0.25 * self.t_c()) * self.pylon_length * self.pylon_chord
        return s_wet_pylon

    def aspect_ratio(self):
        ar_pylon = self.pylon_length ** 2 / self.area()
        return ar_pylon

    """ ------------------------------ Ratios used for normalization --------------------------------------------- """
    def area_ratio(self):
        ar_pylon = self.area() / self.ref_area
        return ar_pylon

    def area_ratio_wet(self):
        ar_w_pylon = self.area_wetted() / self.ref_area
        return ar_w_pylon

    def velocity_ratio(self):
        v_ratio = self.inflow_velocity() ** 2 / self.v_inf ** 2
        return v_ratio

    def chord_ratio(self):
        c_ratio = self.pylon_chord / self.ref_chord
        return c_ratio

    """ ------------------------------ Coefficients for force calculations --------------------------------------- """
    def cl(self):
        """ coefficient taken from Xfoil"""
        cl_polar = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(self.inflow_angle()[0]))
        cl_pylon = float(cl_polar[0]) * np.cos(np.radians(self.cant_angle))

        cl_pylon_norm = cl_pylon * self.area_ratio() * self.velocity_ratio()
        return cl_pylon, cl_pylon_norm

    def cy(self):
        """ coefficient taken from Xfoil"""
        cy_polar = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(self.inflow_angle()[0]))
        cy_pylon = float(cy_polar[0]) * np.sin(np.radians(self.cant_angle))

        cy_pylon_norm = cy_pylon * self.area_ratio() * self.velocity_ratio()
        return cy_pylon, cy_pylon_norm

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)

        ftc = (1 + 1.2 * self.t_c() + 60 * self.t_c() ** 4)

        cd0_val = cf * fm * ftc * self.area_ratio_wet()
        return cd0_val

    def cdi(self):
        cdi_pylon = self.cl()[0] ** 2 / (0.95 * np.pi * self.aspect_ratio())

        cdi_pylon_norm = self.cl()[0] ** 2 / (0.95 * np.pi * self.aspect_ratio())
        return cdi_pylon, cdi_pylon_norm

    def cd(self):
        """ coefficient taken from Xfoil """
        cd_polar = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(self.inflow_angle()[0]))
        cd_value = float(cd_polar[1])

        cd_value_norm = cd_value * self.area_ratio() * self.velocity_ratio()
        return cd_value, cd_value_norm

    def cm(self):
        """ coefficient taken from Xfoil """
        cm_polar = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(self.inflow_angle()[0]))
        cm_pylon_norm = float(cm_polar[2]) * self.velocity_ratio() * self.area_ratio() * self.chord_ratio()
        cm_pylon = float(cm_polar[2])
        return cm_pylon, cm_pylon_norm

    def cn(self):
        alpha = self.inflow_angle()[1]

        cn_pylon = self.cl()[0] * np.cos(alpha) + self.cd()[0] * np.sin(alpha)

        cn_norm = cn_pylon * self.area_ratio() * self.velocity_ratio()
        return cn_pylon, cn_norm

    def ct(self):
        alpha = self.inflow_angle()[1]

        ct_pylon = self.cl()[0] * np.sin(alpha) - self.cd()[0] * np.cos(alpha)

        ct_norm = ct_pylon * self.area_ratio() * self.velocity_ratio()
        return ct_pylon, ct_norm

    """ ---------------------------------- Interference effects----------------------------------------------------- """
    def cd_interference(self):
        cd_pylon_duct = (drag_interference(self.t_c(), "t-junction")
                         * self.area_ratio() * self.velocity_ratio())  # interference pylon duct

        cd_pylon_fuse = (drag_interference(self.t_c(), "plane")
                         * self.area_ratio() * self.velocity_ratio())  # interference pylon fuselage

        cd_int_pylon = cd_pylon_duct + cd_pylon_fuse
        return cd_int_pylon

    """ ---------------------------------- Forces --------------------------------------------------------- """
    def lift_force(self):
        lift_pylon = self.cl()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return lift_pylon

    def drag_force(self):
        drag_pylon = self.cd()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return drag_pylon

    def side_force(self):
        side_pylon = self.cy()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return side_pylon

    def moment_force(self):
        moment_pylon = self.cm()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area() * self.pylon_chord
        return moment_pylon

    """ ------------------------------- Weight definition of the pylon ------------------------------------------- """
    @staticmethod
    def weight():
        """ Pylon weight is calculated in empennage_assembly_PE.py """
        m_pylon = None
        return m_pylon


"""      Test section      """
"""
if __name__ == "__main__":

    a = np.linspace(0, 15, 31)
    cl = []
    #a_ref = np.linspace(0, 15, 16)
    cl_ref = []
    cl_the = []
    cd = []
    cd_ref = []
    cd_the = []
    cd0 = []
    cm = []

    for i in range(len(a)):
        pylon = Pylon(pylon_length=config.pylon_length,
                      pylon_chord=config.pylon_chord,
                      pylon_profile=config.pylon_airfoil,
                      power_condition="on",
                      cant_angle=config.cant_angle,
                      alpha=a[i],
                      ref_area=ref.s_w,
                      v_inf=128,
                      m_supported=10000,
                      ref_chord=2.2345)
        polar = airfoil_polar(f"pylon0012.txt", float(a[i]))
        cd_val = float(polar[1])
        cl_val = float(polar[0])
        cm_val = float(polar[2])

        al = np.radians(a[i])
        cl_theory = np.pi * 2 * al
        cl_the.append(cl_theory)
        cd_the.append(pylon.cd0())

        cl.append(pylon.cl())
        cd.append(pylon.cd())
        cd_ref.append(cd_val)
        cl_ref.append(cl_val)
        cd0.append(pylon.cd0() * pylon.ref_area - 0.002)
        cm.append(pylon.cm()[0])

    plt.figure('CL - alpha')
    plt.plot(a, cl, label=r'Model', color="tab:blue")
    plt.plot(a, cl_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(a, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Pylon')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Model', color="tab:blue")
    plt.plot(a, cd_ref, label=r'XFoil', color="tab:green", marker='o')
    #plt.plot(a, cd_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'Pylon')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD')
    plt.plot(cd, cl, label=r'Model', color="tab:blue")
    plt.plot(cd_ref, cl_ref, label=r'XFoil', color="tab:green", marker='o')
    #plt.plot(cd_the, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Pylon')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - CL^2')
    plt.plot([cl ** 2 for cl in cl], cd, label=r'Model', color="tab:blue")
    plt.plot([cl_ref ** 2 for cl_ref in cl_ref], cd_ref, label='XFoil', color="tab:green", marker='o')
    plt.plot([cl ** 2 for cl in cl], cd0, label=r'$C_{D0}$', color="tab:blue", linestyle="--")
    plt.xlabel(r'$C_{L}^{2}$ [-]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'Pylon')
    #plt.xlim([0, 0.5])
    plt.ylim([0, 0.07])
    plt.legend()
    plt.grid(True)

    plt.figure('CM - alpha')
    plt.plot(a, cm, label=r'Model', color="tab:blue")
    #plt.plot(a_ref, cm_ref, label=r'Experimental', color="tab:green", marker='o')  # for proper comparison change AR of the calculated values
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'Pylon')
    plt.legend()
    plt.grid(True)


    plt.show()


    print(f"inflow vel: {pylon.inflow_velocity()}")
    print(f"inflow ang: {pylon.inflow_angle()}")
    print(f"area: {pylon.area()}, wet area: {pylon.area_wetted()}")
    print(f"cd0: {pylon.cd0()}")
    print(f"cd interference: {pylon.cd_interference():.5f}")
    print(f"cd prime: {pylon.cd_prime():.5f}")
    print(f"cl: {pylon.cl():.5f}")
    print(f"cl prime: {pylon.cl_prime():.5f}")
    print(f"weight per: {pylon.weight()}") """
