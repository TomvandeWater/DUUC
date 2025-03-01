import data.atr_reference
from aircraft.conventional_empennage.components_properties.horizontal_tail import HorizontalTail
from aircraft.conventional_empennage.components_properties.propeller import Propeller
from aircraft.conventional_empennage.components_properties.nacelle import Nacelle
from aircraft.conventional_empennage.components_properties.vertical_tail import VerticalTail
import flow_conditions
import data.atr_reference as ref
from analysis_modules.aerodynamic import tail_volume


class ConventionalEmpennage:
    def __init__(self, ht_span: float, ht_chord: float, ht_profile: str, ht_taper: float,
                 ht_sweep: float, ht_croot: float, vt_span: float, vt_chord: float, vt_profile: str,
                 vt_taper: float, vt_sweep: float, vt_croot: float, tail_type: str, rpm: float,
                 alpha: float, power_condition: str, n_blades: float, prop_diameter: float,
                 hub_diameter: float, prop_airfoil: str, prop_sweep: float, prop_pitch: float,
                 c_root: float, c_tip: float, v_inf: float, nacelle_length: float,
                 nacelle_diameter: float, reynolds: float, mach: float):
        super().__init__()
        self.ht_span = ht_span
        self.ht_chord = ht_chord
        self.ht_profile = ht_profile
        self.ht_taper = ht_taper
        self.ht_sweep = ht_sweep
        self.ht_croot = ht_croot
        self.vt_span = vt_span
        self.vt_chord = vt_chord
        self.vt_profile = vt_profile
        self.vt_taper = vt_taper
        self.vt_sweep = vt_sweep
        self.vt_croot = vt_croot
        self.tail_type = tail_type
        self.rpm = rpm
        self.alpha = alpha
        self.power_condition = power_condition
        self.n_blades = n_blades
        self.prop_diameter = prop_diameter
        self.hub_diameter = hub_diameter
        self.prop_airfoil = prop_airfoil
        self.prop_sweep = prop_sweep
        self.prop_pitch = prop_pitch
        self.c_root = c_root
        self.c_tip = c_tip
        self.v_inf = v_inf
        self.nacelle_length = nacelle_length
        self.nacelle_diameter = nacelle_diameter
        self.reynolds = reynolds
        self.mach = mach

        self.area_ref = data.atr_reference.s_w

        self.ht_tail = HorizontalTail(self.ht_span, self.ht_chord, self.ht_profile,
                                      self.ht_taper, self.ht_sweep, self.ht_croot,
                                      self.alpha, self.v_inf, self.area_ref, self.mach,
                                      self.reynolds)

        self.vt_tail = VerticalTail(self.vt_span, self.vt_chord, self.vt_profile,
                                    self.vt_taper, self.vt_sweep, self.vt_croot,
                                    self.tail_type, self.alpha, self.v_inf, self.area_ref,
                                    self.reynolds, self.mach)

        self.propeller = Propeller(self.rpm, self.alpha, self.power_condition, self.n_blades,
                                   self.prop_diameter, self.hub_diameter, self.prop_airfoil,
                                   self.prop_sweep, self.prop_pitch, self.c_root, self.c_tip,
                                   self.v_inf, self.area_ref, self.reynolds)

        self.nacelle = Nacelle(self.nacelle_length, self.nacelle_diameter, self.v_inf, self.alpha,
                               self.area_ref, self.prop_airfoil, self.n_blades, self.mach,
                               self.reynolds)

    """ ----------------------- vectors for plots ---------------------------------------------------------------- """
    def cd0_vector(self):
        cd0_nacelle = self.nacelle.cd0()
        cd0_ht = self.ht_tail.cd0()
        cd0_vt = self.vt_tail.cd0()
        return [cd0_nacelle, cd0_ht, cd0_vt]

    def cd_interference_vector(self):
        cd_interference_ht = self.ht_tail.cd_interference()
        cd_interference_vt = self.vt_tail.cd_interference()
        cd_interference_nac = self.nacelle.cd_interference() * 2
        return [cd_interference_nac, cd_interference_ht, cd_interference_vt]

    """ ------------------------------- coefficient sums ------------------------------------------------------- """
    def cd_interference(self):
        cd_int_conv = self.ht_tail.cd_interference() + self.nacelle.cd_interference()
        return cd_int_conv

    def cd_prime(self):
        cd_prime_conv = (self.ht_tail.cd_prime() + self.vt_tail.cd_prime() +
                         2 * self.propeller.cd_prime() + 2 * self.nacelle.cd_prime()
                         + self.cd_interference())
        return cd_prime_conv

    def cl_prime(self):
        cl_prime_conv = (self.ht_tail.cl_prime() +
                         2 * self.propeller.cl_prime() + 2 * self.nacelle.cl_prime())
        return cl_prime_conv

    def ht_volume_coefficient(self):
        tv = tail_volume(self.ht_tail.area(), ref.lever_h, self.area_ref, ref.c_root_w)
        return tv

    def vt_volume_coefficient(self):
        tv = tail_volume(self.vt_tail.area(), ref.lever_h, self.area_ref, ref.b_w)
        return tv

    """ ------------------------------------------ Forces ------------------------------------------------------- """
    def thrust(self):
        thrust_conv = 2 * self.propeller.thrust()
        return thrust_conv

    def drag(self):
        drag_conv = self.cd_prime() * self.v_inf ** 2 * self.area_ref * 0.5 * flow_conditions.rho
        return drag_conv

    """ --------------------------------------- Weights ---------------------------------------------------------- """
    def weight_emp(self):
        w_emp = self.ht_tail.weight() + self.vt_tail.weight()
        return w_emp

    def weight_prop(self):
        w_prop = 2 * self.propeller.weight_engine()
        return w_prop

    def weight_nac(self):
        w_nac = 2 * self.nacelle.weight()
        return w_nac

    @staticmethod
    def weight_cv():
        """ function is for complete weight control group"""
        ksc = 0.64
        wto = ref.MTOW * 2.20462

        w_cv = ksc * wto ** 0.75
        return w_cv / 2.20462
