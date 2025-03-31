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
    def __init__(self, pylon_length: float, pylon_chord: float, pylon_profile: str,
                 power_condition: str, cant_angle: float, alpha: float, ref_area: float,
                 v_inf: float, m_supported: float, ref_chord, altitude: float):
        super().__init__()
        self.pylon_length = pylon_length
        self.pylon_chord = pylon_chord
        self.pylon_airfoil = pylon_profile
        self.cant_angle = cant_angle
        self.alpha = alpha
        self.pc = power_condition
        self.ref_area = ref_area
        self.v_inf = v_inf
        self.m_supported = m_supported
        self.ref_chord = ref_chord
        self.altitude = altitude

    """ The inflow speed for the pylon is affected by the outside of the duct
    and also affected by the downwash of the wing. The downwash of the wing and the duct affect
    the inlfow angle of on the pylon"""
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
        """ inflow angle already determined in general propulsive empennage class"""
        cant = np.radians(self.cant_angle)

        inflow_angle = self.alpha * np.cos(cant)
        return inflow_angle

    def reynolds_number(self):
        re_pylon = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.pylon_chord)
        return re_pylon

    """ ------------------------- geometric properties ----------------------------------------------------------- """

    def area(self):
        area_pylon = self.pylon_length * self.pylon_chord
        return area_pylon

    def t_c(self):
        num_list = [int(digit) for digit in self.pylon_airfoil]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def area_wetted(self):
        s_wet_pylon = 2 * (1 + 0.5 * self.t_c()) * self.pylon_length * self.pylon_chord
        return s_wet_pylon

    def aspect_ratio(self):
        ar_pylon = self.pylon_length ** 2 / self.area()
        return ar_pylon

    """ ------------------------------ Coefficients for force calculations --------------------------------------- """

    def cl_da(self):
        cl0 = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(0.0))
        cl0_val = float(cl0[0])
        cl5 = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(10.0))
        cl5_val = float(cl5[0])
        cl_da_pylon = (cl5_val - cl0_val) / 10
        return cl_da_pylon

    def cl(self):
        """ linear section modelled based on Xfoil"""
        cl_pylon = self.cl_da() * self.inflow_angle()
        return cl_pylon

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(0.44)
        area_ratio = self.area_wetted() / self.ref_area

        ftc = (1 + 1.2 * self.t_c() + 60 * self.t_c() ** 4)

        cd0_val = cf * fm * ftc * area_ratio
        return cd0_val

    def cdi(self):
        cdi_pylon = self.cl() ** 2 / (0.95 * np.pi * self.aspect_ratio())
        return cdi_pylon

    def cd(self):
        cd_polar = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(self.inflow_angle()))
        cd_value = float(cd_polar[1]) #+ self.cd0()

        #cd_value = self.cd0() + self.cdi()
        return cd_value

    def cd_interference(self):
        norm_area = (self.t_c() * self.pylon_chord ** 2) / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_pylon_duct = (drag_interference(self.t_c(), "t-junction")
                         * norm_area * norm_speed)  # interference pylon duct

        cd_pylon_fuse = (drag_interference(self.t_c(), "plane")
                         * norm_area * norm_speed)  # interference pylon fuselage

        cd_int_pylon = cd_pylon_duct + cd_pylon_fuse
        return cd_int_pylon

    def coefficient_norm(self):
        norm_area = (self.t_c() * self.pylon_chord ** 2) / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_norm = self.cd() * norm_speed * norm_area
        cl_norm = self.cl() * norm_speed * norm_area
        return cl_norm, cd_norm

    """ ---------------------------------- output primes --------------------------------------------------------- """
    def cl_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cl_pylon = self.cl() * norm_area * norm_speed
        return cl_pylon

    def cd_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_pylon = self.cdi() * norm_area * norm_speed + self.cd0()
        return cd_pylon

    def cm(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_chord = self.pylon_chord / self.ref_chord

        cm_polar = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(self.inflow_angle()))
        cm_pylon_norm = float(cm_polar[2]) * norm_speed * norm_area * norm_chord
        cm_pylon = float(cm_polar[2])
        return cm_pylon, cm_pylon_norm

    def cn(self):
        n_area = self.area() / self.ref_area
        n_velocity = self.inflow_velocity() ** 2 / self.v_inf ** 2

        alpha = np.radians(self.inflow_angle())

        cn_pylon = self.cl() * np.cos(alpha) + self.cd() * np.sin(alpha)

        cn_norm = cn_pylon * n_area * n_velocity
        return cn_pylon, cn_norm

    def ct(self):
        n_area = self.area() / self.ref_area
        n_velocity = self.inflow_velocity() ** 2 / self.v_inf ** 2

        alpha = np.radians(self.inflow_angle())

        ct_pylon = self.cl() * np.sin(alpha) - self.cd() * np.cos(alpha)

        ct_norm = ct_pylon * n_area * n_velocity
        return ct_pylon, ct_norm

    """ ------------------------------- Weight definition of the pylon ------------------------------------------- """
    def weight(self):
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
