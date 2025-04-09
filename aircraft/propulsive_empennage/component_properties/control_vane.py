import numpy as np
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa
from data.read_data import airfoil_polar
from analysis_modules.factors import k_control, skin_friction, mach_correction, oswald
import data.atr_reference as ref
import config


class ControlVane:
    """ Reference for each control vane is the attachment to the nacelle the values calculated
    are for 1 control vane. """
    def __init__(self, geometry, conditions, reference, power_condition: str, va_inlet: float,
                 d_exit: float, deflection: float, v_after_prop: float, a_after_prop: float):
        super().__init__()
        self.cv_span = geometry[0]
        self.cv_chord = geometry[1]
        self.cv_profile = geometry[2]

        self.pc = power_condition

        self.ref_area = reference[0]
        self.ref_chord = reference[1]

        self.v_inf = conditions[0]
        self.alpha = conditions[1]
        self.altitude = conditions[2]
        self.mach = conditions[3]
        self.density = air_density_isa(self.altitude)[0]

        self.v_after_prop = v_after_prop
        self.a_after_prop = a_after_prop
        self.va_inlet = va_inlet
        self.d_exit = d_exit
        self.deflection_angle = deflection

    """ ------------------------------ Determine inflow properties ------------------------------------------------ """
    def inflow_velocity(self):
        area_vane = np.pi / 4 * self.d_exit ** 2
        v_vane = self.va_inlet / area_vane

        if self.pc == "off":
            v_cv = v_vane
            return v_cv
        else:
            v_cv = np.cos(np.radians(self.a_after_prop / 2)) * self.v_after_prop
            return v_cv

    def inflow_angle(self):
        if self.pc == "off":
            alpha_power_off = 0
            deflection = self.deflection_angle
            inflow = alpha_power_off + deflection

            angle_rad = np.radians(inflow)
            return inflow, angle_rad
        else:
            inflow = (self.a_after_prop / 2) + self.deflection_angle
            angle_rad = np.radians(inflow)
            return inflow, angle_rad

    def reynolds_number(self):
        re_cv = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.cv_chord)
        return re_cv

    """ -------------------------------- geometric properties --------------------------------------------------- """
    def area(self):
        area_control_vane = self.cv_span * self.cv_chord
        return area_control_vane

    def aspect_ratio(self):
        aspect_ratio_control = self.cv_span ** 2 / self.area()
        return aspect_ratio_control

    def t_c(self):
        num_list = [int(digit) for digit in self.cv_profile]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def area_wetted(self):
        wet_vane = 2 * (1 + 0.25 * self.t_c()) * self.cv_span * self.cv_chord
        return wet_vane

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
        c_ratio = self.cv_chord / self.ref_chord
        return c_ratio
    """ ---------------------------------- determine coefficients based on theory -------------------------------- """
    def cl(self):
        coeff = airfoil_polar(f"vcv{self.cv_profile}.txt", self.inflow_angle()[0])
        cl_vane = float(coeff[0])

        cl_norm = cl_vane * self.area_ratio() * self.velocity_ratio()
        return cl_vane, cl_norm

    def cd0(self):
        cf = skin_friction(self.reynolds_number(), "t")
        fm = mach_correction(self.mach)
        ftc = 1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4

        coeff = airfoil_polar(f"vcv{self.cv_profile}.txt", float(0.0))
        cdmin = float(coeff[1])

        cd0_vane = cf * fm * ftc * (cdmin / 0.004) ** 0.4
        cd0_vane_ar = cf * fm * ftc * self.area_ratio_wet() * (cdmin / 0.004) ** 0.4
        return cd0_vane, cd0_vane_ar

    def cdi(self):
        e = oswald(self.aspect_ratio(), 0)
        cdi_vane = self.cl()[0] ** 2 / (np.pi * self.aspect_ratio() * e)
        cdi_vane_norm = (self.cl()[0] ** 2 / (np.pi * self.aspect_ratio() * e) * self.velocity_ratio()
                         * self.area_ratio())
        return cdi_vane, cdi_vane_norm

    def cd(self):
        coeff = airfoil_polar(f"vcv{self.cv_profile}.txt", self.inflow_angle()[0])
        cd_vane = float(coeff[1])
        cd_norm = cd_vane * self.velocity_ratio() * self.area_ratio()
        return cd_vane, cd_norm

    def cm(self):
        """ defined at xac = 0.25c for symmetrical airfoil"""
        coeff = airfoil_polar(f"vcv{self.cv_profile}.txt", self.inflow_angle()[0])
        cm_vane = float(coeff[2])

        cm_norm = cm_vane * self.area_ratio() * self.velocity_ratio() * self.chord_ratio()
        return cm_vane, cm_norm

    def cn(self):
        alpha = self.inflow_angle()[1]

        cn_pylon = self.cl()[0] * np.cos(alpha) + self.cd()[0] * np.sin(alpha)

        cn_norm = cn_pylon * self.area_ratio() * self.velocity_ratio()
        return cn_pylon, cn_norm

    def ct(self):
        alpha = self.inflow_angle()[1]

        ct_pylon = self.cl()[0] * np.sin(alpha) - self.cd()[0] * np.cos(alpha)

        ct_norm = ct_pylon * self.area_ratio() * self.velocity_ratio()
        return ct_pylon, ct_norm

    """ ------------------------------ INTERFERENCE EFFECTS  ------------------------------------------------------- """
    def cd_interference(self):
        # assumed constant chord and thickness for the control vanes
        cd_cv_nac = (drag_interference(self.t_c(), "plane")
                     * self.area_ratio() * self.velocity_ratio())

        cd_cv_duct = (drag_interference(self.t_c(), "t-junction")
                      * self.area_ratio() * self.velocity_ratio())

        cd_interference_vcv = cd_cv_nac + cd_cv_duct
        return cd_interference_vcv

    """ ------------------------------ FORCES  --------------------------------------------------------------------- """
    def lift_force(self):
        lift_control = self.cl()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return lift_control

    def drag_force(self):
        drag_control = self.cd()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area()
        return drag_control

    def moment_force(self):
        moment_control = self.cm()[0] * 0.5 * self.density * self.inflow_velocity() ** 2 * self.area() * self.cv_chord
        return moment_control

    """ -------------------------------- WEIGHT -------------------------------------------------------------------- """
    @staticmethod
    def weight():
        """ function is for complete weight control group"""
        ksc = 0.64
        wto = ref.MTOW * 2.20462

        # weight division 0.75 wing - 0.25 tail
        w_cv = (ksc * wto ** (2/3)) / 32  # assume equally spread over the 4 control surfaces and both PE's
        return w_cv / 2.20462


""" Test section"""
"""
if __name__ == "__main__":
    control = ControlVane(cv_span=config.control_vane_length,
                          cv_chord=config.control_vane_chord,
                          cv_profile=config.control_vanes_airfoil,
                          power_condition="on",
                          va_inlet=128 * (np.pi * (config.duct_diameter / 2)),
                          d_exit=config.d_exit,
                          u1=130,
                          ref_area=config.duct_diameter * config.duct_chord,
                          deflection=0,
                          v_inf=128,
                          alpha=0,
                          mach=0.576, u_mom=10, a_after_prop=50)

    print(f"inflow vel: {control.inflow_velocity()}")
    print(f"inflow ang: {control.inflow_angle()}")
    print(f"cd: {control.cd():.5f}")
    print(f"cd0: {control.cd0():.5f}")
    print(f"cdi: {control.cdi():.5f}")
    print(f"cd interference: {control.cd_interference():.5f}")
    print(f"cd prime: {control.cd_prime():.5f}")
    print(f"cl: {control.cl():.5f}")
    print(f"cl prime: {control.cl_prime():.5f}")
    print(f"weight: {control.weight()}")
    print(f"wet area:{control.wet_area()}") """
