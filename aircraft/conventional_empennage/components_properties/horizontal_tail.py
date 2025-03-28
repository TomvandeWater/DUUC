import numpy as np
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa
from analysis_modules.factors import skin_friction, mach_correction, oswald
from data.read_data import airfoil_polar
import data.atr_reference as ref
import matplotlib.pyplot as plt


class HorizontalTail:
    """ Horizontal tail class for the conventional empennage """
    def __init__(self, ht_span: float, ht_chord: float, ht_profile: str,
                 ht_taper: float, ht_sweep: float, ht_croot: float, alpha: float, v_inf: float,
                 area_ref: float, mach: float, ar_wing: float, cl_wing: float, cla_wing: float,
                 altitude: float):
        super().__init__()
        self.ht_span = ht_span
        self.ht_chord = ht_chord
        self.ht_profile = ht_profile
        self.ht_taper = ht_taper
        self.ht_sweep = ht_sweep
        self.ht_croot = ht_croot
        self.alpha = alpha
        self.v_inf = v_inf
        self.area_ref = area_ref
        self.mach = mach
        self.ar_wing = ar_wing
        self.cl_wing = cl_wing
        self.cla_wing = cla_wing
        self.altitude = altitude

    """ ---------------------------- Calculate inflow properties ------------------------------------------------- """
    def inflow_velocity(self):
        de = (2 * self.cl_wing) / (np.pi * self.ar_wing)
        inflow_ht = self.v_inf * (1 - (de ** 2) / 2)
        return inflow_ht

    def inflow_angle(self):
        e0 = (2 * self.cl_wing) / (np.pi * self.ar_wing)
        de_da = 2 * self.cla_wing / (np.pi * self.ar_wing)

        eta = e0 + de_da * self.alpha

        inflow_angle = self.alpha - ref.alpha_install_wing - eta + ref.installation_angle
        # inflow_angle = self.alpha -> for plotting validation
        return inflow_angle

    def reynolds_number(self):
        re_htail = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.ht_chord)
        return re_htail

    """" ----------------------- Calculate geometric properties of the horizontal tail ---------------------------- """
    def tip_chord(self):
        c_tip = self.ht_croot * self.ht_taper
        return c_tip

    @staticmethod
    def area():
        s_ht = ref.s_ht
        return s_ht

    def t_c(self):
        num_list = [int(digit) for digit in self.ht_profile]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def wet_area(self):
        wet_ht = 2 * (1 + 0.5 * self.t_c()) * self.ht_span * self.ht_chord
        return wet_ht

    def aspect_ratio(self):
        aspect_ratio_ht = self.ht_span ** 2 / self.area()
        return aspect_ratio_ht

    def x_ht_tail(self):
        """ Taking the leading edge of the chord profile as a reference"""
        sweep_rad = np.radians(self.ht_sweep)
        x_tip_tail = 0.5 * self.ht_span * np.sin(sweep_rad)

        return x_tip_tail

    def oswald(self):
        f = 0.005 * (1 + 1.5 * (ref.tr_h - 0.6) ** 2)
        mc = 1 + 0.12 * self.mach ** 6
        c1 = (0.142 + f * self.aspect_ratio() * (10 * self.t_c()) ** 0.33) / (np.cos(np.radians(ref.phi_qc_h)) ** 2)
        c2 = (0.1 * (3 * 0 + 1)) / (4 + self.aspect_ratio()) ** 0.8

        e_ht = 1 / (mc * (1 + c1 + c2))

        e_ht = oswald(self.aspect_ratio(), 0)  # used improved oswald factor calculation
        return e_ht

    """ ------------------------------------ determine coefficients -------------------------------------------  """
    def cd_interference(self):
        norm_area = (self.t_c() * self.ht_chord ** 2) / self.area()
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_int_ht = (drag_interference(self.t_c(), 't-junction') * norm_speed * norm_area)

        return cd_int_ht

    def cl_al(self):
        tan2 = np.tan(np.radians(self.ht_sweep)) ** 2
        an = 2 * np.pi * self.aspect_ratio()
        bn = 2 + np.sqrt(4 + self.aspect_ratio() ** 2 * (1 + tan2))
        cl_al_tail = an / bn

        cl_al_tail = (2 * np.pi * self.aspect_ratio()) / (2 + np.sqrt(4 + self.aspect_ratio() ** 2))
        return cl_al_tail

    def cl(self):
        cl_polar = airfoil_polar(f"ht{self.ht_profile}.txt", 0)
        cl_ht0 = cl_polar[0]

        cl_ht = cl_ht0 + self.cl_al() * np.radians(self.inflow_angle())
        return cl_ht

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)
        sweep_corr = 1.34 * self.mach ** 0.18 * (np.cos(np.radians(self.ht_sweep)) ** 0.28)
        ftc = (1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4) * sweep_corr
        norm_area = self.area() / self.area_ref

        cd0_ht = cf * fm * ftc * 1.04 * norm_area   # 1.04 for conventional empennage talking about interference
        # cd0_ht = cf * fm * ftc * 1.04
        return cd0_ht

    def cdi(self):
        e = self.oswald()
        cdi_ht = (self.cl() ** 2) / (np.pi * self.aspect_ratio() * e)
        return cdi_ht

    def cd(self):
        norm_area = self.area() / self.area_ref
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_ht = self.cd0() * norm_speed + self.cdi() * norm_speed * norm_area
        cd_ht = self.cd0() + self.cdi()
        return cd_ht

    @staticmethod
    def cm():
        """ as this contribution is small to the whole aircraft we assume it zero"""
        cm_ht = 0
        return cm_ht

    """ ------------------------ Output primes ----------------------------------------------------------------- """
    def cd_prime(self):
        norm_area = self.area() / self.area_ref
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_prime_ht = self.cd() * norm_speed * norm_area
        return cd_prime_ht

    def cl_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.area() / self.area_ref

        cl_prime_vt = self.cl() * norm_speed * norm_area
        return cl_prime_vt

    def ct(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.area() / self.area_ref

        alpha = np.radians(self.inflow_angle())

        ct = self.cl() * np.sin(alpha) - self.cd() * np.cos(alpha)
        ct_norm = ct * norm_area * norm_speed

        return ct, ct_norm

    def cn(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.area() / self.area_ref

        alpha = np.radians(self.inflow_angle())

        cn = self.cl() * np.cos(alpha) + self.cd() * np.sin(alpha)
        cn_norm = cn * norm_area * norm_speed

        return cn, cn_norm

    """" ----------------------------------- Calculate weight --------------------------------------------------- """
    def weight(self):
        kh = 1.1
        sh = self.area()
        vd = ref.v_dive
        sweep = np.radians(ref.phi_hc_h)

        m_hor = kh * sh * (62 * (sh ** 0.2 * vd) / (1000 * np.sqrt(np.cos(sweep))) - 2.5)
        return m_hor


""" Test section"""

if __name__ == "__main__":

    a = np.linspace(0, 15, 31)
    cl = []
    #a_ref = np.linspace(0, 15, 16)
    cl_ref = []
    cl_the = []
    cd = []
    cd_ref = []
    cd_the = []

    for i in range(len(a)):

        hor = HorizontalTail(ht_span=ref.b_h,
                             ht_chord=ref.c_root_h,
                             ht_profile=ref.airfoil_ht,
                             ht_taper=ref.tr_h,
                             ht_sweep=ref.phi_qc_h,
                             ht_croot=ref.c_root_h,
                             alpha=a[i],
                             v_inf=128,
                             area_ref=ref.s_w,
                             mach=0.443,
                             ar_wing=12,
                             cl_wing=1.44,
                             cla_wing=5.89
                             )

        polar = airfoil_polar(f"ht0009.txt", float(a[i]))
        cd_val = float(polar[1])
        cl_val = float(polar[0])

        al = np.radians(a[i])
        cl_theory = 5.7 * al  # from experimental data
        cl_the.append(cl_theory)
        cd_the.append(hor.cd0() + (cl_theory ** 2) / (np.pi * hor.aspect_ratio() * 0.6))

        cl.append(hor.cl())
        cd.append(hor.cd())
        cd_ref.append(cd_val)
        cl_ref.append(cl_val)

    plt.figure('CL - alpha')
    plt.plot(a, cl, label=r'Model', color="tab:orange")
    plt.plot(a, cl_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(a, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Horizontal tail')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Model', color="tab:orange")
    plt.plot(a, cd_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(a, cd_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'Horizontal tail')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD')
    plt.plot(cd, cl, label=r'Model', color="tab:orange")
    plt.plot(cd_ref, cl_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(cd_the, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Horizontal tail')
    plt.legend()
    plt.grid(True)

    [cl_val ** 2 for cl_val in cl]

    plt.figure('CD - CL^2')
    plt.plot([cl ** 2 for cl in cl], cd, label=r'Model', color="tab:orange")
    plt.plot([cl_ref ** 2 for cl_ref in cl_ref], cd_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot([cl_the ** 2 for cl_the in cl_the], cd_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$C_{L}^2$ [-]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'Horizontal tail')
    plt.xlim(0, 1)
    plt.ylim(0, 0.15)
    plt.legend()
    plt.grid(True)

    plt.show()


    print(f"inflow vel: {hor.inflow_velocity()}")
    print(f"inflow ang: {hor.inflow_angle()}")
    print(f"cd0: {hor.cd0()}, cd: {hor.cd():.5f}")
    print(f"cd prime: {hor.cd_prime():.5f}")
    print(f"cl: {hor.cl():.5f}")
    print(f"cl prime: {hor.cl_prime():.5f}")
    print(f"weight: {hor.weight()}")
    print(f"area: {hor.area()}")
    print(f"span: {hor.ht_span}")
    print(f"sweep: {hor.ht_sweep}")
    print(f"ar: {hor.aspect_ratio()}")
    print(f"area: {hor.area()}, wet_area: {hor.wet_area()}")
