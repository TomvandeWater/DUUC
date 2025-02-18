import numpy as np
from analysis_modules.factors import oswald, skin_friction, mach_correction
from data.read_data import airfoil_polar


class VerticalTail:
    """ Vertical Tail class for conventional empennage"""
    def __init__(self, vt_span: float, vt_chord: float, vt_profile: str,
                 vt_taper: float, vt_sweep: float, vt_croot: float, tail_type: str,
                 alpha: float, v_inf: float, area_ref: float, reynolds: float, mach: float):
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
        self.reynolds = reynolds

    def inflow_velocity(self):
        inflow_vt = self.v_inf
        return inflow_vt

    @staticmethod
    def inflow_angle():
        inflow_angle = 0
        return inflow_angle

    """" Calculate geometric properties of the vertical tail"""

    def tip_chord(self):
        c_tip = self.vt_croot * self.vt_taper
        return c_tip

    def area(self):
        """ Area from the vertical stabilizer (based on trapezoid area)"""

        s_vt = 0.5 * self.vt_span * 0.5 * (self.vt_croot + self.vt_chord)
        print(f'Area of the vertical stabilizer = {2* s_vt} [m^2]')
        return s_vt

    def t_c(self):
        num_list = [int(digit) for digit in self.vt_profile]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def wet_area(self):
        wet_vt = 2 * (1 + 0.5 * self.t_c()) * self.vt_span * self.vt_chord
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

    def cl(self):
        cl_vt = airfoil_polar(f"vt{self.vt_profile}.txt", self.inflow_angle())
        return cl_vt

    def cdi(self):
        e = oswald(self.aspect_ratio(), self.vt_sweep)
        cdi_vt = self.cl() ** 2 / (np.pi * self.aspect_ratio() * e)
        return cdi_vt

    def cd0(self):
        cf = skin_friction(self.reynolds, "t")
        fm = mach_correction(self.mach)
        norm_area = self.wet_area() / self.area_ref
        ftc = 1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4

        coeff = airfoil_polar(f"vt{self.vt_profile}.txt", float(0.0))
        cdmin = float(coeff[1] + coeff[2])
        cd0_vt = cf * fm * ftc * norm_area * (cdmin / 0.004) ** 0.4
        return cd0_vt

    def cd(self):
        cd_vt = self.cd0() + self.cdi()
        return cd_vt

    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.area() / self.area_ref

        alpha = np.radians(self.inflow_angle())

        cd_cd = self.cd() * np.cos(alpha) * norm_speed * norm_area
        cd_cl = self.cl() * np.sin(alpha) * norm_speed * norm_area

        cd_prime_vt = cd_cl + cd_cd
        return cd_prime_vt

    def c_beta_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.area() / self.area_ref

        alpha = np.radians(self.inflow_angle())

        cl_cl = self.cl() * np.cos(alpha) * norm_speed * norm_area
        cl_cd = self.cd() * np.sin(alpha) * norm_speed * norm_area

        cl_prime_vt = cl_cl + cl_cd
        return cl_prime_vt

