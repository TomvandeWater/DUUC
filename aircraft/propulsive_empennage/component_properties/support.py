import numpy as np
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa
from data.read_data import airfoil_polar
from analysis_modules.factors import skin_friction, mach_correction


class SupportStrut:
    def __init__(self, geometry, conditions, reference, power_condition: str, v_after_prop: float, va_inlet: float,
                 prop_diameter: float, a_after_prop: float):
        super().__init__()
        self.support_length = geometry[0]
        self.support_chord = geometry[1]
        self.support_profile = geometry[2]
        self.cant_angle = geometry[3]

        self.pc = power_condition

        self.v_inf = conditions[0]
        self.alpha = conditions[1]
        self.altitude = conditions[2]
        self.mach = conditions[3]
        self.density = air_density_isa(self.altitude)[0]

        self.ref_area = reference[0]
        self.ref_chord = reference[1]

        self.va_inlet = va_inlet
        self.prop_diameter = prop_diameter
        self.a_after_prop = a_after_prop
        self.v_after_prop = v_after_prop

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
        if self.pc == "off":
            cant = np.radians(self.cant_angle)
            inflow_support = - cant * np.cos(0)
            angle_rad = np.radians(inflow_support)
            return inflow_support, angle_rad
        else:
            cant = np.radians(self.cant_angle)

            inflow_support = self.a_after_prop + - cant * np.cos(0)
            angle_rad = np.radians(inflow_support)
            return inflow_support, angle_rad

    def reynolds_number(self):
        re_sup = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.support_chord)
        return re_sup

    """  ------------------------------------- geometric properties -------------------------------------------- """
    def area(self):
        """ flat surface area of the support """
        s_support = self.support_chord * self.support_length
        return s_support

    def area_wetted(self):
        """ untapered, unswept wing definition by Torenbeek """
        s_wet_support = 2 * (1 + 0.25 * self.t_c()) * self.support_length * self.support_chord
        return s_wet_support

    def t_c(self):
        """ assume NACA 4-series airfoil """
        num_list = [int(digit) for digit in self.support_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def aspect_ratio(self):
        ar_support = self.support_length ** 2 / self.area()
        return ar_support

    """ ------------------------------ Ratios used for normalization --------------------------------------------- """
    def area_ratio(self):
        ar_support = self.area() / self.ref_area
        return ar_support

    def area_ratio_wet(self):
        """ Area ratio used in zero lift drag """
        ar_w_support = self.area_wetted() / self.ref_area
        return ar_w_support

    def velocity_ratio(self):
        v_ratio = self.inflow_velocity() ** 2 / self.v_inf ** 2
        return v_ratio

    def chord_ratio(self):
        c_ratio = self.support_chord / self.ref_chord
        return c_ratio

    """ -------------------------------------- coefficients for support forces ------------------------------------ """
    def cl(self):
        """ coefficient taken from Xfoil"""
        cl_polar = airfoil_polar(f"support{self.support_profile}.txt", float(self.inflow_angle()[0]))
        cl_support = float(cl_polar[0]) * np.cos(np.radians(self.cant_angle))

        cl_support_norm = cl_support * self.velocity_ratio() * self.area_ratio()
        return cl_support, cl_support_norm

    def cy(self):
        """ coefficient taken from Xfoil"""
        cy_polar = airfoil_polar(f"support{self.support_profile}.txt", float(self.inflow_angle()[0]))
        cy_support = float(cy_polar[0]) * np.sin(np.radians(self.cant_angle))

        cy_support_norm = cy_support * self.velocity_ratio() * self.area_ratio()
        return cy_support, cy_support_norm


    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)

        ftc = (1 + 1.2 * self.t_c() + 60 * self.t_c() ** 4)

        cd0_support = cf * fm * ftc
        cd0_support_ratio = cf * fm * ftc * self.area_ratio_wet()
        return cd0_support, cd0_support_ratio

    def cdi(self):
        cdi_support = self.cl()[0] ** 2 / (0.95 * np.pi * self.aspect_ratio())

        cdi_norm = self.cl()[0] ** 2 / (0.95 * np.pi * self.aspect_ratio()) * self.area_ratio() * self.velocity_ratio()
        return cdi_support, cdi_norm

    def cd(self):
        """ coefficient taken from Xfoil"""
        cd_polar = airfoil_polar(f"support{self.support_profile}.txt", float(self.inflow_angle()[0]))
        cd_val = float(cd_polar[1])

        cd_norm = cd_val * self.area_ratio() * self.velocity_ratio()
        return cd_val, cd_norm

    def cm(self):
        """ coefficient taken from Xfoil"""
        cm_strut = airfoil_polar(f"support{self.support_profile}.txt", float(self.inflow_angle()[0]))
        cm_strut = float(cm_strut[2])
        cm_norm = cm_strut * self.area_ratio() * self.velocity_ratio() * self.chord_ratio()
        return cm_strut, cm_norm

    def cn(self):
        alpha = self.inflow_angle()[1]

        cn_support = self.cl()[0] * np.cos(alpha) + self.cd()[0] * np.sin(alpha)

        cn_norm = cn_support * self.area_ratio() * self.velocity_ratio()
        return cn_support, cn_norm

    def ct(self):
        alpha =  self.inflow_angle()[1]

        ct_support = self.cl()[0] * np.sin(alpha) - self.cd()[0] * np.cos(alpha)

        ct_norm = ct_support * self.area_ratio() * self.velocity_ratio()

        return ct_support, ct_norm

    """ ---------------------------------- Interference effects----------------------------------------------------- """
    def cd_interference(self):
        cd_support_nacelle = (2 * drag_interference(self.t_c(), "plane")
                              * self.area_ratio() * self.velocity_ratio())  # multiplied by 2 for 2 nac-supp inter

        cd_support_duct = (2 * drag_interference(self.t_c(), "t-junction")
                           * self.area_ratio() * self.velocity_ratio())  # multiplied by 2 for 2 supp-duct inter

        cd_int_support = cd_support_nacelle + cd_support_duct
        return cd_int_support

    """ ------------------------------- Forces ------------------------------------------------- """
    def lift_force(self):
        lift_support = self.cl()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return lift_support

    def drag_force(self):
        drag_support = self.cd()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return drag_support

    def moment_force(self):
        moment_support = (self.cm()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
                          * self.support_chord)
        return moment_support

    def side_force(self):
        side_support = self.cy()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return side_support

    """ ----------------------------------- Weight estimation ------------------------------------------------- """
    @staticmethod
    def weight():
        """ Weight is calculated in the empennage_assembly_PE.py file"""
        m_support = None
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
                               ref_area=ref.s_w,
                               v_inf=128,
                               a_after_prop=a[i],
                               m_supported=1000,
                               prop_diameter=3.6,
                               va_inlet=250,
                               ref_chord=2.2345)

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

    print(f"inflow vel: {support.inflow_velocity()}")
    print(f"inflow ang: {support.inflow_angle()}")
    print(f"cd0: {support.cd0()}")
    print(f"cd interference: {support.cd_interference():.5f}")
    print(f"cd prime: {support.cd_prime():.5f}")
    print(f"cl: {support.cl():.5f}")
    print(f"cl prime: {support.cl_prime():.5f}")
    print(f"wet area: {support.area_wetted()}") 
    """
