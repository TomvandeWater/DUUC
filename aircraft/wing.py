import numpy as np
import data.atr_reference as ref
from analysis_modules.factors import skin_friction, mach_correction
from data.read_data import airfoil_polar
import matplotlib.pyplot as plt


class Wing:
    def __init__(self, b_wing: float, sweep_wing: float, wing_airfoil: str, taper_ratio_wing: float, c_root_wing: float,
                 alpha: float, velocity: float, mach: float, reynolds_number: float, wing_type: str):
        super().__init__()
        self.b_wing = b_wing
        self.sweep_wing = sweep_wing
        self.wing_airfoil = wing_airfoil
        self.tr_wing = taper_ratio_wing
        self.c_root = c_root_wing
        self.alpha = alpha
        self.velocity = velocity
        self.mach = mach
        self.re = reynolds_number
        self.wing_type = wing_type

    """ -------------------------------- Define inflow properties ------------------------------------------------ """
    def inflow_velocity(self):
        v_wing = self.velocity
        return v_wing

    def inflow_angle(self):
        """ inflow angle in degrees"""
        a_wing = self.alpha
        return a_wing

    """ ------------------------------------ Define geometric properties ---------------------------------------- """
    @staticmethod
    def area():
        area_wing = ref.s_w
        return area_wing

    def aspect_ratio(self):
        ar_wing = self.b_wing ** 2 / self.area()
        return ar_wing

    def t_c(self):
        """ assume wing airfoil is a NACA 5 series airfoil"""
        num_list = [int(digit) for digit in self.wing_airfoil]
        thickness = num_list[3] * 10 + num_list[4]  # naca thickness of profile
        t_c_wing = thickness / 100  # returns value in percentage of normalized chord

        x_camb = (num_list[1] * 10 + num_list[2]) / 100

        return t_c_wing, x_camb

    def area_wetted(self):
        s_wet = 2 * (1 + 0.5 * (self.t_c()[0])) * self.b_wing * self.c_root
        return s_wet

    def oswald(self):
        e_wing = (1.78 * (1 - 0.0045 * self.aspect_ratio() ** 0.68)) / np.sqrt(self.aspect_ratio())
        if self.wing_type == "conventional":
            n_eng = 2
        else:
            n_eng = 0

        f = 0.005 * (1 + 1.5 * (self.tr_wing - 0.6) ** 2)
        mc = 1 + 0.12 * self.mach ** 6
        c1 = (0.142 + f * self.aspect_ratio() * (10 * self.t_c()[0]) ** 0.33)/(np.cos(np.radians(ref.phi_qc_w)) ** 2)
        c2 = (0.1 * (3 * n_eng + 1)) / (4 + self.aspect_ratio()) ** 0.8

        e_wing = 1 / (mc * (1 + c1 + c2))
        return e_wing

    """ ------------------------------------ determine coefficients -------------------------------------------- """
    def cd0(self):
        fm = mach_correction(self.mach)
        cf = skin_friction(self.re, "t")
        f_wing = ((1 + (0.6 / 0.15) * self.t_c()[0] + 100 * self.t_c()[0] ** 4)
                  * (1.34 * self.mach ** 0.18 * (np.cos(np.radians(ref.phi_qc_w))) ** 0.28))

        cd0_wing = fm * cf * f_wing
        return cd0_wing

    def cdi(self):
        cdi_wing = self.cl() ** 2 / (np.pi * self.aspect_ratio() * self.oswald())
        return cdi_wing

    def cd(self):
        cd_wing = self.cd0() + self.cdi()
        return cd_wing

    def cl0(self):
        cl0_w = airfoil_polar(f"wing{self.wing_airfoil}.txt", float(ref.alpha_install_wing))
        cl0_wing = cl0_w[0]
        return cl0_wing

    def cl_al(self):
        cl_al_wing = (np.pi * 2 * self.aspect_ratio()) / (2 + np.sqrt(4 + self.aspect_ratio() ** 2
                                                                      * (1 + np.tan(np.radians(self.sweep_wing)) ** 2
                                                                         - self.mach ** 2)))
        return cl_al_wing

    def cl(self):
        if self.wing_type == "conventional":
            angle = self.inflow_angle() - 1
        else:
            angle = self.inflow_angle()

        cl_wing = self.cl0() + self.cl_al() * np.radians(angle)
        return cl_wing

    def cm0(self):
        cm0_wing = -0.1 * (self.t_c()[1] / self.c_root)

        sweep = np.radians(self.sweep_wing)

        sweep_corr = np.cos(sweep) ** 2

        cm_wing = cm0_wing * sweep_corr
        return cm_wing

    def cm_a(self):
        cma_wing = - 0.05 / np.sqrt(self.aspect_ratio())
        return cma_wing

    def cm(self):
        cm_wing = self.cm0() + self.cm_a() * self.alpha
        return cm_wing

    """ ---------------------------- determine output primes ----------------------------------------------------- """
    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.velocity ** 2

        alpha = np.radians(self.inflow_angle())

        cd_cd = self.cd() * np.cos(alpha) * norm_speed

        cd_cl = self.cl() * np.sin(alpha) * norm_speed

        cd_prime_wing = cd_cd + cd_cl
        return cd_prime_wing

    def cl_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.velocity ** 2

        alpha = np.radians(self.inflow_angle())

        cl_cl = self.cl() * np.cos(alpha) * norm_speed

        cl_cd = self.cd() * np.sin(alpha) * norm_speed

        cl_prime_wing = cl_cl + cl_cd
        return cl_prime_wing

    """ ------------------------------------ determine weight of the wing --------------------------------------- """
    def weight(self):
        """ based on Torenbeek class II weight estimation"""
        w_zf = ref.ZFW_max * 9.81  # convert kg to N
        kw = 6.67e-3
        n_limit = 2.5
        n_ult = 1.5 * n_limit
        bs = ref.b_w / np.cos(np.radians(ref.phi_qc_w))
        tr = self.t_c()[0]

        sqrt_b = np.sqrt(ref.b_w / bs)
        ratio2 = ((bs / tr) / (w_zf / ref.s_w))

        if self.wing_type == "conventional":

            w_wing_r = kw * bs ** 0.75 * (1 + sqrt_b) * n_ult ** 0.55 * ratio2 ** 0.30
            w_wing = w_wing_r * w_zf  # returns wing weight in Newtons
            return ((1 - 0.03) * w_wing) / 9.81

        if self.wing_type == "DUUC":
            f_no_engine = 0.85  # reduced wing weight due to no engines on the wing between 0.75-0.85

            w_wing_r = kw * bs ** 0.75 * (1 + sqrt_b) * n_ult ** 0.55 * ratio2 ** 0.30 * f_no_engine
            w_wing = w_wing_r * w_zf  # returns wing weight in Newtons

            return (1.02 * w_wing) / 9.81

        else:
            print("Wrong input for wing type")
            return None


""" Test section """

if __name__ == "__main__":

    a = np.linspace(0, 15, 31)
    axfoil = np.linspace(0, 11, 23)
    cl = []
    cl_duc = []
    cd_duc = []
    #a_ref = np.linspace(0, 15, 16)
    cl_ref = []
    cl_the = []
    cd = []
    cd_ref = []
    cd_the = []

    for i in range(len(axfoil)):
        polar = airfoil_polar(f"wing43013.txt", float(a[i] + ref.alpha_install_wing))
        cd_val = float(polar[1] + polar[2])
        cl_val = float(polar[0])
        cd_ref.append(cd_val)
        cl_ref.append(cl_val)

    for i in range(len(a)):

        wing = Wing(b_wing=ref.b_w,
                    sweep_wing=ref.phi_qc_w,
                    wing_airfoil=ref.wing_airfoil,
                    taper_ratio_wing=ref.tr_w,
                    c_root_wing=ref.c_root_w,
                    alpha=a[i],
                    velocity=128,
                    mach=0.443,
                    reynolds_number=8422274,
                    wing_type="conventional")

        duc = Wing(b_wing=ref.b_w,
                   sweep_wing=ref.phi_qc_w,
                   wing_airfoil=ref.wing_airfoil,
                   taper_ratio_wing=ref.tr_w,
                   c_root_wing=ref.c_root_w,
                   alpha=a[i],
                   velocity=128,
                   mach=0.443,
                   reynolds_number=8422274,
                   wing_type="DUUC")

        al = np.radians(a[i] + ref.alpha_install_wing)
        cl_theory = np.pi * 2 * al
        cl_the.append(cl_theory)
        cd_the.append(wing.cd0() + cl_theory ** 2 / (np.pi * wing.aspect_ratio() * 0.77))

        cl.append(wing.cl())
        cd.append(wing.cd())
        cl_duc.append(duc.cl())
        cd_duc.append(duc.cd())


    plt.figure('CL - alpha')
    plt.plot(a, cl, label=r'Model - ATR72-600', color="tab:orange")
    plt.plot(a, cl_duc, label=r'Model - DUUC', color="tab:blue")
    plt.plot(axfoil, cl_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(a, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Wing')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha')
    plt.plot(a, cd, label=r'Model - ATR72-600', color="tab:orange")
    plt.plot(a, cd_duc, label=r'Model - DUUC', color="tab:blue")
    plt.plot(axfoil, cd_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(a, cd_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'Wing')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD')
    plt.plot(cd, cl, label=r'Model  - ATR72-600', color="tab:orange")
    plt.plot(cd_duc, cl_duc, label=r'Model - DUUC', color="tab:blue")
    plt.plot(cd_ref, cl_ref, label=r'XFoil', color="tab:green", marker='o')
    plt.plot(cd_the, cl_the, label='Theory', color="tab:green", linestyle="--")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'Wing')
    plt.legend()
    plt.grid(True)
    plt.show()


    print(f"inflow vel: {wing.inflow_velocity()}")
    print(f"inflow ang: {wing.inflow_angle()}, oswald: {wing.oswald()}")
    print(f"area: {wing.area()}, wetted area: {wing.area_wetted()}")
    print(f"aspect ratio: {wing.aspect_ratio()}, t_c: {wing.t_c()}")
    print(f"cd0: {wing.cd0()}, cdi: {wing.cdi()}, cd: {wing.cd()}, cdprime: {wing.cd_prime()}")
    print(f"cl0: {wing.cl0()}, cl_al: {wing.cl_al()}, cl: {wing.cl()}, cl_prime: {wing.cl_prime()}")
    print(f"cm0: {wing.cm0()}, cma: {wing.cm_a()}, cm: {wing.cm()}")
    print(f"weight: {wing.weight()}")
