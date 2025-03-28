import numpy as np
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.factors import mach_correction, skin_friction
from analysis_modules.ISA import air_density_isa


class Nacelle:
    """ Nacelle class for conventional configuration"""
    def __init__(self, nacelle_length: float, nacelle_diameter: float, v_inf: float, alpha: float,
                 area_ref: float, prop_airfoil: str, n_blades: float, mach: float, altitude: float):
        super().__init__()
        self.nacelle_length = nacelle_length
        self.nacelle_diameter = nacelle_diameter
        self.v_inf = v_inf
        self.alpha = alpha
        self.area_ref = area_ref
        self.prop_airfoil = prop_airfoil
        self.n_blades = n_blades
        self.mach = mach
        self.altitude = altitude

    """ ---------------------------- Calculate inflow properties ------------------------------------------------- """
    def inflow_velocity(self):
        inflow_nac = self.v_inf
        return inflow_nac

    def inflow_angle(self):
        alfa = self.alpha
        return alfa

    def reynolds_number(self):
        re_nac = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.nacelle_length)
        return re_nac

    """" --------------------------- Calculate geometric properties ----------------------------------------------- """

    def area_wetted(self):
        """ only the side area is assumed no closing sides"""
        area_cylinder = np.pi * self.nacelle_diameter * self.nacelle_length

        """ half a sphere is used to close of the rear end"""
        area_rear = 0.5 * np.pi * self.nacelle_diameter ** 2

        area_wet_nac = area_rear + area_cylinder
        return area_wet_nac

    def area(self):
        area_nac = self.nacelle_length * self.nacelle_diameter
        return area_nac

    def t_c(self):
        """ defined on comparable propeller blade"""
        t_c_prop_root = 0.10
        t_c_prop_tip = 0.08
        t_c_av = (t_c_prop_tip + t_c_prop_root) / 2
        return t_c_av

    """" ----------------------- Calculate geometric properties of the horizontal tail ---------------------------- """
    def cd_interference(self):
        norm_area = self.area() / self.area_ref
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_int_nac = (drag_interference(self.t_c(), 'plane') * norm_speed * norm_area
                      * self.n_blades)
        return cd_int_nac

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)
        l_d = self.nacelle_length / self.nacelle_diameter
        # f_nac = 1 + 60 / l_d ** 3 + 0.0025 * l_d
        f_nac = 1 + 0.35 / (self.nacelle_length / self.nacelle_diameter)

        cd0_nacelle = cf * fm * f_nac * 1.5 * (self.area() / self.area_ref)
        return cd0_nacelle

    def cd(self):
        cd_nac = self.cd0() * (self.area_ref / self.area())
        return cd_nac

    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2
        norm_area = self.area() / self.area_ref

        cd_nac = self.cd0() * norm_speed * norm_area
        return cd_nac

    @staticmethod
    def cl():
        """assume the nacelle does not produce lift"""
        cl_prime_nac = 0
        return cl_prime_nac

    """ ------------------------ Output primes ----------------------------------------------------------------- """
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
    @staticmethod
    def weight():
        """ based on Torenbeek Class II weight estimation"""
        # calculated for take off conditions
        p = 4102000
        eta = 0.64545
        v = 54.12

        m_nacelle = 0.0458 * (p * eta) / (v * 9.81)
        return m_nacelle


""" Test section"""
"""
if __name__ == "__main__":
    hor = Nacelle(nacelle_length=ref.l_nacelle,
                  nacelle_diameter=ref.d_nacelle,
                  alpha=0,
                  v_inf=128,
                  area_ref=ref.s_w,
                  mach=0.576,
                  reynolds=8422274,
                  n_blades=ref.n_blades,
                  prop_airfoil=ref.propeller_airfoil)

    print(f"inflow vel: {hor.inflow_velocity()}")
    print(f"cd: {hor.cd0():.5f}, cd0: {hor.cd0()}")
    print(f"cd prime: {hor.cd_prime():.5f}")
    print(f"cl: {hor.cl_prime():.5f}")

    print(f"weight: {hor.weight()}")
    print(f"area: {hor.area()}, wet: {hor.area_wetted()}")"""
