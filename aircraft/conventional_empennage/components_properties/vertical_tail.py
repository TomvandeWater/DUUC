import numpy as np
from analysis_modules.factors import oswald, skin_friction, mach_correction
from data.read_data import airfoil_polar
import data.atr_reference as ref
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa


class VerticalTail:
    """ Vertical Tail class for conventional empennage"""
    def __init__(self, vt_span: float, vt_chord: float, vt_profile: str,
                 vt_taper: float, vt_sweep: float, vt_croot: float, tail_type: str,
                 alpha: float, v_inf: float, area_ref: float, mach: float, altitude: float):
        super().__init__()
        self.vt_span = vt_span
        self.vt_chord = vt_chord
        self.vt_profile = vt_profile
        self.tail_type = tail_type
        self.vt_taper = vt_taper
        self.vt_sweep = vt_sweep
        self.vt_croot = vt_croot
        self.v_inf = v_inf
        self.alpha = alpha
        self.area_ref = area_ref
        self.mach = mach
        self.altitude = altitude

    """ ---------------------------- Calculate inflow properties ------------------------------------------------- """
    def inflow_velocity(self):
        inflow_vt = self.v_inf
        return inflow_vt

    @staticmethod
    def inflow_angle():
        inflow_angle = 0
        return inflow_angle

    def reynolds_number(self):
        re_vtail = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.vt_chord)
        return re_vtail

    """" --------------------------- Calculate geometric properties of the vertical tail -------------------------- """

    def tip_chord(self):
        c_tip = self.vt_croot * self.vt_taper
        return c_tip

    @staticmethod
    def area():
        """ Area from the vertical stabilizer (based on trapezoid area)"""
        s_vt = ref.s_vt
        return s_vt

    def t_c(self):
        num_list = [int(digit) for digit in self.vt_profile]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def wet_area(self):
        wet_vt = 2 * self.area() * (1 + 0.25 * self.t_c() * (1 + 0.7 * ref.tr_v) / (1 + ref.tr_v))
        return wet_vt

    def aspect_ratio(self):
        aspect_ratio_vt = self.vt_span ** 2 / self.area()
        return aspect_ratio_vt

    def x_vt_tail(self):
        """ Taking the leading edge of the chord profile as a reference"""
        sweep_rad = np.radians(self.vt_sweep)
        x_tip_tail = 0.5 * self.vt_span * np.sin(sweep_rad)
        return x_tip_tail

    def z_vt_tail(self):
        """ depending on the tail type determine the z-location of the tail"""
        if self.tail_type == 'T':
            z_tail = self.vt_span
            return z_tail

        else:
            z_tail = 0
            return z_tail

    """ -------------------------------- Coefficient calculation -------------------------------------------------- """
    def cd_interference(self):
        norm_area = (self.t_c() * self.vt_chord ** 2) / self.area()
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_int_vt = (drag_interference(self.t_c(), 'plane') * norm_speed
                     * norm_area)
        return cd_int_vt

    def cl(self):
        cl_polar = airfoil_polar(f"vt{self.vt_profile}.txt", self.inflow_angle())
        cl_vt = cl_polar[0]
        return cl_vt

    def cdi(self):
        e = oswald(self.aspect_ratio(), self.vt_sweep)
        cdi_vt = self.cl() ** 2 / (np.pi * self.aspect_ratio() * e)
        return cdi_vt

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)
        norm_area = self.wet_area() / self.area_ref
        sweep_corr = 1.34 * self.mach ** 0.18 * (np.cos(np.radians(self.vt_sweep)) ** 0.28)
        ftc = (1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4) * sweep_corr

        cd0_vt = cf * fm * ftc * norm_area * 1.04
        return cd0_vt

    def cd(self):
        norm_area = self.wet_area() / self.area_ref
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_vt = self.cd0() * norm_speed + self.cdi() * norm_speed * norm_area
        return cd_vt

    """ ------------------------------ calculate output primes ---------------------------------------------------- """
    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        alpha = np.radians(self.inflow_angle())

        cd_cd = self.cd() * np.cos(alpha) * norm_speed
        cd_cl = self.cl() * np.sin(alpha) * norm_speed

        cd_prime_vt = cd_cl + cd_cd
        return cd_prime_vt

    @staticmethod
    def cl_prime():
        """ vertical tail does not produce lift"""
        cl_prime_vt = 0
        return cl_prime_vt

    def c_beta_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.area() / self.area_ref

        alpha = np.radians(self.inflow_angle())

        cl_cl = self.cl() * np.cos(alpha) * norm_speed * norm_area
        cl_cd = self.cd() * np.sin(alpha) * norm_speed * norm_area

        cl_prime_vt = cl_cl + cl_cd
        return cl_prime_vt

    """ ---------------------------------------- determine weight ------------------------------------------------- """
    def weight(self):
        kv = 1
        sv = self.area()
        vd = ref.v_dive

        sweep = np.radians(ref.phi_hc_v)

        m_hor = kv * sv * (62 * (sv ** 0.2 * vd) / (1000 * np.sqrt(np.cos(sweep))) - 2.5)
        return m_hor


""" Test section"""

if __name__ == "__main__":
    hor = VerticalTail(vt_span=ref.b_v,
                       vt_chord=ref.c_root_v,
                       vt_profile=ref.airfoil_vt,
                       vt_taper=ref.tr_v,
                       vt_sweep=ref.phi_qc_v,
                       vt_croot=ref.c_root_v,
                       alpha=0,
                       v_inf=128,
                       area_ref=ref.s_w,
                       mach=0.44, tail_type="t-tail")

    print(f"inflow vel: {hor.inflow_velocity()}")
    print(f"inflow ang: {hor.inflow_angle()}")
    print(f"cd: {hor.cd():.5f}")
    print(f"cd prime: {hor.cd_prime():.5f}")
    print(f"cd0: {hor.cd0()}")
    print(f"cl: {hor.cl():.5f}")

    print(f"weight: {hor.weight()}")
    print(f"area: {hor.area()}")
    print(f"span: {hor.vt_span}")
    print(f"area: {hor.area()}, wet_area: {hor.wet_area()}")
