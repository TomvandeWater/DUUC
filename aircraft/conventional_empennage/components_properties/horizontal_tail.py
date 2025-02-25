import numpy as np
from analysis_modules.aerodynamic import drag_interference
from analysis_modules.factors import skin_friction, mach_correction, oswald
from data.read_data import airfoil_polar
import data.atr_reference as ref


class HorizontalTail:
    """ Horizontal tail class for the conventional empennage """
    def __init__(self, ht_span: float, ht_chord: float, ht_profile: str,
                 ht_taper: float, ht_sweep: float, ht_croot: float, alpha: float, v_inf: float,
                 area_ref: float, mach: float, reynolds: float):
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
        self.reynolds = reynolds

    def inflow_velocity(self):
        de = (2 * ref.cl_wing) / (np.pi * ref.ar_w)
        inflow_ht = self.v_inf * (1 - (de ** 2) / 2)
        return inflow_ht

    def inflow_angle(self):
        de = (2 * ref.cl_wing) / (np.pi * ref.ar_w)
        inflow_angle = self.alpha - de
        return inflow_angle

    """" Calculate geometric properties of the horizontal tail"""

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

    def cd_int(self):
        norm_area = (self.t_c() * self.ht_chord ** 2) / self.area()
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_int_ht = (drag_interference(self.t_c(), 't-junction') * norm_speed
                     * norm_area * 2)
        return cd_int_ht

    def cl(self):
        cl_polar = airfoil_polar(f"ht{self.ht_profile}.txt", self.inflow_angle())
        cl_ht = cl_polar[0]
        return cl_ht

    def cd0(self):
        cf = skin_friction(self.reynolds, "t")
        fm = mach_correction(self.mach)
        norm_area = self.wet_area() / self.area_ref
        ftc = 1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4

        coeff = airfoil_polar(f"ht{self.ht_profile}.txt", float(0.0))
        cdmin = float(coeff[1] + coeff[2])
        cd0_ht = cf * fm * ftc * norm_area * (cdmin / 0.004) ** 0.4
        return cd0_ht

    def cdi(self):
        e = oswald(self.aspect_ratio(), 0)
        cdi_ht = self.cl() ** 2 / (np.pi * self.aspect_ratio() * e)
        return cdi_ht

    def cd(self):
        cd_ht = self.cd0() + self.cdi()
        return cd_ht

    def cd_prime(self):
        norm_area = self.area() / self.area_ref
        norm_speed = self.inflow_angle() ** 2 / self.v_inf ** 2

        alpha = np.radians(self.inflow_angle())

        cd_cd = self.cd() * np.cos(alpha) * norm_speed * norm_area
        cd_cl = self.cl() * np.sin(alpha) * norm_speed * norm_area

        cd_prime_ht = cd_cd + cd_cl
        return cd_prime_ht

    def cl_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.area() / self.area_ref

        alpha = np.radians(self.inflow_angle())

        cl_cl = self.cl() * np.cos(alpha) * norm_speed * norm_area
        cl_cd = self.cd() * np.sin(alpha) * norm_speed * norm_area

        cl_prime_vt = cl_cl + cl_cd
        return cl_prime_vt

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
    hor = HorizontalTail(ht_span=ref.b_h,
                         ht_chord=ref.c_root_h,
                         ht_profile=ref.airfoil_ht,
                         ht_taper=ref.tr_h,
                         ht_sweep=ref.phi_qc_h,
                         ht_croot=ref.c_root_h,
                         alpha=0,
                         v_inf=128,
                         area_ref=ref.s_w,
                         mach=0.576,
                         reynolds=8422274)

    print(f"inflow vel: {hor.inflow_velocity()}")
    print(f"inflow ang: {hor.inflow_angle()}")
    print(f"cd: {hor.cd():.5f}")
    print(f"cd prime: {hor.cd_prime():.5f}")
    print(f"cl: {hor.cl():.5f}")
    print(f"cl prime: {hor.cl_prime():.5f}")
    print(f"weight: {hor.weight()}")
    print(f"area: {hor.area()}")
    print(f"span: {hor.ht_span}")
    print(f"sweep: {hor.ht_sweep}") """
