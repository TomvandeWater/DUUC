import data.atr_reference
from aircraft.conventional_empennage.components_properties.horizontal_tail import HorizontalTail
from aircraft.conventional_empennage.components_properties.propeller import Propeller
from aircraft.conventional_empennage.components_properties.nacelle import Nacelle
from aircraft.conventional_empennage.components_properties.vertical_tail import VerticalTail
from analysis_modules.ISA import air_density_isa
import data.atr_reference as ref
from analysis_modules.aerodynamic import tail_volume
import matplotlib.pyplot as plt
import numpy as np
import config


class ConventionalEmpennage:
    def __init__(self, ht_span: float, ht_chord: float, ht_profile: str, ht_taper: float,
                 ht_sweep: float, ht_croot: float, vt_span: float, vt_chord: float, vt_profile: str,
                 vt_taper: float, vt_sweep: float, vt_croot: float, tail_type: str, rpm: float,
                 alpha: float, power_condition: str, n_blades: float, prop_diameter: float,
                 hub_diameter: float, prop_airfoil: str, prop_sweep: float, prop_pitch: float,
                 c_root: float, c_tip: float, v_inf: float, nacelle_length: float, nacelle_diameter: float,
                 mach: float, ar_wing: float, cl_wing: float, cla_wing: float, bem_input,
                 altitude: float):
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
        self.mach = mach
        self.ar_wing = ar_wing
        self.cl_wing = cl_wing
        self.cla_wing = cla_wing
        self.bem_input = bem_input
        self.altitude = altitude
        self.area_ref = data.atr_reference.s_w

        """ Inflow angle for htail is corrected for the wing inflow by ar_wing, cl_wing and cla_wing"""
        self.ht_tail = HorizontalTail(ht_span=self.ht_span, ht_chord=self.ht_chord, ht_profile=self.ht_profile,
                                      ht_taper=self.ht_taper, ht_sweep=self.ht_sweep, ht_croot=self.ht_croot,
                                      alpha=self.alpha,  area_ref=self.area_ref, mach=self.mach,
                                      ar_wing=self.ar_wing, cl_wing=self.cl_wing, cla_wing=self.cla_wing,
                                      altitude=self.altitude, v_inf=self.v_inf)

        """ Vtail outputs cl_a in the form of a sideslip coefficient """
        self.vt_tail = VerticalTail(vt_span=self.vt_span, vt_chord=self.vt_chord, vt_profile=self.vt_profile,
                                    vt_taper=self.vt_taper, vt_sweep=self.vt_sweep, vt_croot=self.vt_croot,
                                    tail_type=self.tail_type, alpha=self.alpha, v_inf=self.v_inf,
                                    area_ref=self.area_ref, mach=self.mach, altitude=self.altitude)

        self.propeller = Propeller(self.rpm, self.alpha, self.power_condition, self.n_blades,
                                   self.prop_diameter, self.hub_diameter, self.prop_airfoil,
                                   self.prop_sweep, self.prop_pitch, self.c_root, self.c_tip,
                                   self.v_inf, self.area_ref, altitude=self.altitude)

        self.nacelle = Nacelle(self.nacelle_length, self.nacelle_diameter, self.v_inf, self.alpha,
                               self.area_ref, self.prop_airfoil, self.n_blades, self.mach,
                               altitude=self.altitude)
    """ -------------------------------- geometric properties ------------------------------------------------------ """
    def area_wet(self):
        s_wet_htail = self.ht_tail.wet_area()
        s_wet_vtail = self.vt_tail.wet_area()

        s_wet_tot = s_wet_htail + s_wet_vtail
        return s_wet_tot

    """ ----------------------- vectors for plots ------------------------------------------------------------------ """
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
        """ normalized components to wing area"""
        cd_prime_conv = (self.ht_tail.cd_prime() + self.vt_tail.cd_prime() +
                         2 * self.propeller.cd_prime() + 2 * self.nacelle.cd_prime()
                         + self.cd_interference())
        return cd_prime_conv

    def cl_prime(self):
        """ normalized components to wing area"""
        cl_prime_conv = self.ht_tail.cl_prime() + 2 * self.propeller.cl_prime()
        return cl_prime_conv

    def cl_a(self):
        """ assume htail is only producing significant cl_a - > per radian"""
        cl_a_emp = self.ht_tail.cl_al()
        return cl_a_emp

    def cd_sum(self):
        """ regular cd-values """
        cd_htail = self.ht_tail.cd()
        cd_vtail = self.vt_tail.cd()
        return cd_htail + cd_vtail

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
        drag_conv = self.cd_prime() * self.v_inf ** 2 * self.area_ref * 0.5 * air_density_isa(self.altitude)[0]
        return drag_conv

    def lift(self):
        lift_conv = self.cl_prime() * self.v_inf ** 2 * self.area_ref * 0.5 * air_density_isa(self.altitude)[0]
        return lift_conv

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


""" Test section"""
"""
if __name__ == "__main__":

    a = np.linspace(0, 20, 41)
    cl = []
    cd = []

    # cl_duct = []
    # cd_duct = []
    cl_duct = [0.03897488655623439, 0.0493637629729071, 0.05975263938957981, 0.07014151580625252, 0.08053039222292523, 0.09091926863959793, 0.10130814505627064, 0.11169702147294337, 0.12208589788961606, 0.1324747743062888, 0.14286365072296148, 0.15325252713963422, 0.1636414035563069, 0.1740302799729796, 0.18441915638965234, 0.19480803280632503, 0.20519690922299771, 0.21558578563967046, 0.2259746620563432, 0.2363635384730159, 0.24675241488968858, 0.25714129130636126, 0.267530167723034, 0.2779190441397067, 0.2883079205563794, 0.29869679697305207, 0.3090856733897248, 0.3194745498063975, 0.3298634262230703, 0.340252302639743, 0.3506411790564156, 0.36103005547308836, 0.37141893188976105, 0.3818078083064338, 0.3921966847231065, 0.40258556113977917, 0.41297443755645197, 0.42336331397312466, 0.43375219038979734, 0.44414106680647, 0.4545299432231427]
    cd_duct = [0.0009664262563692861, 0.0010054675006860897, 0.0011220801560485128, 0.001316581877125317, 0.0015887973244081682, 0.0019388049806360224, 0.0023667424712764985, 0.002872281951062508, 0.0034557881617842128, 0.004116865783551538, 0.004855800554822508, 0.005672395788912263, 0.006566832535342762, 0.007539166848709743, 0.008589158191802098, 0.0097172315844463, 0.010922876388136125, 0.012207044755665145, 0.01356897378805938, 0.01500973481434311, 0.016529208733492128, 0.018128318085430924, 0.01981251753144595, 0.021574288388506602, 0.02341363065661287, 0.025330544335764756, 0.02734682489699482, 0.029455414927819272, 0.0316510286719346, 0.03394767492151673, 0.03632189258214446, 0.038775813811733735, 0.04130735489216147, 0.04391751494640508, 0.046605766153527205, 0.049339420610586, 0.05209320570270551, 0.05492456220587062, 0.057833490120081335, 0.06081998944533767, 0.06388406018163967]

    for i in range(len(a)):
        emp = ConventionalEmpennage(ht_span=ref.b_h,
                                    ht_chord=ref.c_root_h,
                                    ht_profile=ref.airfoil_ht,
                                    ht_taper=ref.tr_h,
                                    ht_sweep=ref.phi_qc_h,
                                    ht_croot=ref.c_root_h,
                                    vt_span=ref.b_v,
                                    vt_chord=ref.c_root_v,
                                    vt_profile=ref.airfoil_vt,
                                    vt_taper=ref.tr_v,
                                    vt_sweep=ref.phi_qc_v,
                                    vt_croot=ref.c_root_v,
                                    tail_type="t-tail",
                                    rpm=config.rpm,
                                    alpha=a[i],
                                    power_condition=config.power_condition,
                                    n_blades=config.n_blades,
                                    prop_diameter=config.duct_diameter * 0.95,
                                    hub_diameter=config.hub_diameter,
                                    prop_airfoil=config.prop_airfoil,
                                    prop_sweep=config.propeller_sweep,
                                    prop_pitch=config.propeller_pitch,
                                    c_root=config.c_root,
                                    c_tip=config.c_tip,
                                    nacelle_length=ref.l_nacelle,
                                    nacelle_diameter=ref.d_nacelle,
                                    v_inf=128,
                                    mach=0.44,
                                    ar_wing=12,
                                    cl_wing=1.44,
                                    cla_wing=5.89,
                                    altitude=7000)

        cl.append(emp.cl_prime())
        cd.append(emp.cd_sum())

    plt.figure('CL - alpha - Empennage comparison')
    plt.plot(a, cl, label=r'T-tail ATR72-600', color="tab:blue")
    plt.plot(a, cl_duct, label=r'Propulsive empennage', color="tab:orange")
    plt.xlim([0, 20])
    # plt.ylim([0, 0.5])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $\alpha$')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha - Empennage comparison')
    plt.plot(a, cd, label=r'T-tail ATR72-600', color="tab:blue")
    plt.plot(a, cd_duct, label=r'Propulsive empennage', color="tab:orange")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $\alpha$')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD - Empennage comparison')
    plt.plot(cd, cl, label=r'T-tail ATR72-600', color="tab:blue")
    plt.plot(cd_duct, cl_duct, label=r'Propulsive empennage', color="tab:orange")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    # plt.xlim([0, 0.04])
    # plt.ylim([-0.1, 0.3])
    plt.title(r'$C_{L}$ vs. $C_{D}$')
    plt.legend()
    plt.grid(True)
    plt.show() """
