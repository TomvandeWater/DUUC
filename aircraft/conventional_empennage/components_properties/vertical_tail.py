import numpy as np
from analysis_modules.factors import oswald, skin_friction, mach_correction
from data.read_data import airfoil_polar
import data.atr_reference as ref
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa


class VerticalTail:
    """ Vertical Tail class for conventional empennage"""
    def __init__(self, geometry, conditions, reference, tail_type: str):
        super().__init__()
        self.vt_span = geometry[0]
        self.vt_chord = geometry[1]
        self.vt_profile = geometry[2]
        self.vt_sweep = geometry[3]
        self.vt_taper = geometry[4]
        self.vt_croot = geometry[5]

        self.tail_type = tail_type

        self.v_inf = conditions[0]
        self.alpha = conditions[1]
        self.altitude = conditions[2]
        self.mach = conditions[3]
        self.density = air_density_isa(self.altitude)[0]

        self.ref_area = reference[0]
        self.ref_chord = reference[1]

    """ ---------------------------- Calculate inflow properties ------------------------------------------------- """
    def inflow_velocity(self):
        inflow_vt = self.v_inf
        return inflow_vt

    @staticmethod
    def inflow_angle():
        inflow_angle = 0
        angle_rad = np.radians(inflow_angle)
        return inflow_angle, angle_rad

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

    def area_wetted(self):
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
        c_ratio = self.vt_chord / self.ref_chord
        return c_ratio

    """ -------------------------------- Coefficient calculation -------------------------------------------------- """
    def cl(self):
        cl_polar = airfoil_polar(f"vt{self.vt_profile}.txt", self.inflow_angle()[0])
        cl_vt = cl_polar[0]

        cl_norm = cl_vt * self.area_ratio() * self.velocity_ratio()
        return cl_vt, cl_norm

    def cdi(self):
        e = oswald(self.aspect_ratio(), self.vt_sweep)
        cdi_vt = self.cl()[0] ** 2 / (np.pi * self.aspect_ratio() * e)

        cdi_norm = cdi_vt * self.area_ratio() * self.velocity_ratio()
        return cdi_vt, cdi_norm

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)
        sweep_corr = 1.34 * self.mach ** 0.18 * (np.cos(np.radians(self.vt_sweep)) ** 0.28)
        ftc = (1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4) * sweep_corr

        cd0_vt_ar = cf * fm * ftc * self.area_ratio_wet() * 1.04
        cd0_vt = cf * fm * ftc * 1.04
        return cd0_vt, cd0_vt_ar

    def cd(self):
        cd_vt_norm = self.cd0()[0] * self.area_ratio() + self.cdi()[0] * self.area_ratio() * self.velocity_ratio()
        cd_vt = self.cd0()[1] + self.cdi()[0]
        return cd_vt, cd_vt_norm

    def c_beta(self):
        alpha = self.inflow_angle()[1]

        cbeta_vt = self.cl()[0] * np.cos(alpha) + self.cd()[0] * np.sin(alpha)

        cbeta_norm = cbeta_vt * self.area_ratio() * self.velocity_ratio()
        return cbeta_vt, cbeta_norm

    def ct(self):
        alpha = self.inflow_angle()[1]

        ct_vt = self.cl()[0] * np.sin(alpha) - self.cd()[0] * np.cos(alpha)

        ct_norm = ct_vt * self.area_ratio() * self.velocity_ratio()
        return ct_vt, ct_norm

    def cm(self):
        """ coefficient taken from Xfoil """
        cm_polar = airfoil_polar(f"vt{self.vt_profile}.txt", float(self.inflow_angle()[0]))
        cm_pylon_norm = float(cm_polar[2]) * self.velocity_ratio() * self.area_ratio() * self.chord_ratio()
        cm_pylon = float(cm_polar[2])
        return cm_pylon, cm_pylon_norm

    """ -------------------------------- INTERFERENCE EFFECTS ------------------------------------------------------ """
    def cd_interference(self):
        cd_int_vt = (drag_interference(self.t_c(), 'plane')
                     * self.area_ratio() * self.velocity_ratio())
        return cd_int_vt

    """ ---------------------------------------- FORCES ------------------------------------------------------------ """
    def side_force(self):
        side_force_vt = self.cl()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return side_force_vt

    def drag_force(self):
        drag_vt = self.cd()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return drag_vt

    def moment_force(self):
        moment_vt = self.cm()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area() * self.vt_chord
        return moment_vt

    """ ---------------------------------------------- WEIGHT ------------------------------------------------------ """
    def weight(self):
        kv = 1
        sv = self.area()
        vd = ref.v_dive

        sweep = np.radians(ref.phi_hc_v)

        m_hor = kv * sv * (62 * (sv ** 0.2 * vd) / (1000 * np.sqrt(np.cos(sweep))) - 2.5)
        return m_hor


""" Test section"""
"""
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
    """
