import numpy as np
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa
from analysis_modules.factors import skin_friction, mach_correction, oswald
from data.read_data import airfoil_polar
import data.atr_reference as ref
import matplotlib.pyplot as plt


class HorizontalTail:
    """ Horizontal tail class for the conventional empennage """
    def __init__(self, geometry, conditions, reference, ar_wing: float, cl_wing: float, cla_wing: float):
        super().__init__()
        self.ht_span = geometry[0]
        self.ht_chord = geometry[1]
        self.ht_profile = geometry[2]

        self.ht_sweep = geometry[4]
        self.ht_taper = geometry[3]
        self.ht_croot = geometry[5]

        self.v_inf = conditions[0]
        self.alpha = conditions[1]
        self.altitude = conditions[2]
        self.mach = conditions[3]
        self.density = air_density_isa(self.altitude)[0]

        self.ref_area = reference[0]
        self.ref_chord = reference[1]

        self.ar_wing = ar_wing
        self.cl_wing = cl_wing
        self.cla_wing = cla_wing

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
        angle_rad = np.radians(inflow_angle)
        return inflow_angle, angle_rad

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

    def area_wetted(self):
        wet_ht = 2 * (1 + 0.25 * self.t_c() * ((1 + 0.7 * 0.6) / (1 + 0.6))) * self.area()
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
        e_ht = oswald(self.aspect_ratio(), 0)
        return e_ht

    """ ------------------------------ RATIONS USED FOR NORMALALIZATION ------------------------------------------ """
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
        c_ratio = self.ht_chord / self.ref_chord
        return c_ratio

    """ ------------------------------------ COEFFICIENTS --------------------------------------------------------  """
    def cl_al(self):
        tan2 = np.tan(np.radians(self.ht_sweep)) ** 2

        cl_al_tail = (2 * np.pi * self.aspect_ratio()) / (2 + np.sqrt(4 + self.aspect_ratio() ** 2))
        return cl_al_tail

    def cl(self):
        cl_polar = airfoil_polar(f"ht{self.ht_profile}.txt", 0)
        cl_ht0 = cl_polar[0]

        cl_ht = cl_ht0 + self.cl_al() * self.inflow_angle()[1]
        cl_norm = cl_ht * self.area_ratio() * self.velocity_ratio()
        return cl_ht, cl_norm

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)
        sweep_corr = 1.34 * self.mach ** 0.18 * (np.cos(np.radians(self.ht_sweep)) ** 0.28)
        ftc = (1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4) * sweep_corr

        cd0_ht_ar = cf * fm * ftc * 1.04 * self.area_ratio_wet()
        cd0_ht = cf * fm * ftc * 1.04
        return cd0_ht, cd0_ht_ar

    def cdi(self):
        e = self.oswald()
        cdi_ht = (self.cl()[0] ** 2) / (np.pi * self.aspect_ratio() * e)
        cdi_norm = cdi_ht * self.area_ratio() * self.velocity_ratio()
        return cdi_ht, cdi_norm

    def cd(self):
        cd_ht_norm = self.cd0()[1] * self.velocity_ratio() + self.cdi()[1]
        cd_ht = self.cd0()[1] + self.cdi()[0]
        return cd_ht, cd_ht_norm

    @staticmethod
    def cm():
        """ as this contribution is small to the whole aircraft we assume it zero"""
        cm_ht = 0
        return cm_ht

    def ct(self):
        alpha = self.inflow_angle()[1]

        ct = self.cl()[0] * np.sin(alpha) - self.cd()[0] * np.cos(alpha)
        ct_norm = ct * self.velocity_ratio() * self.area_ratio()
        return ct, ct_norm

    def cn(self):
        alpha = self.inflow_angle()[1]

        cn = self.cl()[0] * np.cos(alpha) + self.cd()[0] * np.sin(alpha)
        cn_norm = cn * self.velocity_ratio() * self.area_ratio()

        return cn, cn_norm

    """ ------------------------------------ INTERFERENCE EFFECTS -------------------------------------------------  """
    def cd_interference(self):
        norm_area = (self.t_c() * self.ht_chord ** 2) / self.area()
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_int_ht = (drag_interference(self.t_c(), 't-junction') * norm_speed * norm_area)

        return cd_int_ht

    """ ------------------------------------- FORCES --------------------------------------------------------------- """
    def lift_force(self):
        lift_htail = self.cl()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return lift_htail

    def drag_force(self):
        drag_htail = self.cd()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return drag_htail

    def moment_force(self):
        moment_htail = self.cm() * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area() * self.ht_chord
        return moment_htail

    """" ------------------------------------- WEIGHT -------------------------------------------------------------- """
    def weight(self):
        kh = 1.1
        sh = self.area()
        vd = ref.v_dive
        sweep = np.radians(ref.phi_hc_h)

        m_hor = kh * sh * (62 * (sh ** 0.2 * vd) / (1000 * np.sqrt(np.cos(sweep))) - 2.5)
        return m_hor


""" Test section"""
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

    for i in range(len(a)):

        hor = HorizontalTail(geometry=[ref.b_h, ref.c_root_h, ref.airfoil_ht, ref.tr_h, ref.phi_qc_h, ref.c_root_h],
                             conditions=[128, a[i], 7000, 0.41], reference=[ref.s_w, 2.45],
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

        cl.append(hor.cl()[0])
        cd.append(hor.cd()[0])
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
    print(f"cd0: {hor.cd0()}, cd: {hor.cd()}")
    print(f"cl: {hor.cl()}")
    print(f"weight: {hor.weight()}")
    print(f"area: {hor.area()}")
    print(f"span: {hor.ht_span}")
    print(f"sweep: {hor.ht_sweep}")
    print(f"ar: {hor.aspect_ratio()}")
    print(f"area: {hor.area()}, wet_area: {hor.area_wetted()}")
    """
