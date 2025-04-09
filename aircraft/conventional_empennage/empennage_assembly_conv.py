import data.atr_reference
from aircraft.conventional_empennage.components_properties.horizontal_tail import HorizontalTail
from aircraft.conventional_empennage.components_properties.propeller import Propeller
from aircraft.conventional_empennage.components_properties.nacelle import Nacelle
from aircraft.conventional_empennage.components_properties.vertical_tail import VerticalTail
import data.atr_reference as ref
from analysis_modules.aerodynamic import tail_volume
import matplotlib.pyplot as plt


class ConventionalEmpennage:
    def __init__(self, ht_geometry, vt_geometry, nacelle_geometry, conditions, reference, power_condition: str,
                 ar_wing: float, cl_wing: float, cla_wing: float, bem_input, rpm: float, geometry_prop):
        super().__init__()
        self.rpm = rpm
        self.power_condition = power_condition

        self.ar_wing = ar_wing
        self.cl_wing = cl_wing
        self.cla_wing = cla_wing
        self.bem_input = bem_input

        self.conditions = conditions
        self.reference = reference
        self.ht_geometry = ht_geometry
        self.vt_geometry = vt_geometry
        self.nacelle_geometry = nacelle_geometry
        self.prop_geometry = geometry_prop

        """ ---------------------------------------- COMPONENT DEFINITIONS ----------------------------------------- """
        self.ht_tail = HorizontalTail(geometry=self.ht_geometry, conditions=self.conditions, reference=self.reference,
                                      ar_wing=self.ar_wing, cl_wing=self.cl_wing, cla_wing=self.cla_wing)

        self.vt_tail = VerticalTail(geometry=self.vt_geometry, conditions=self.conditions, reference=self.reference,
                                    tail_type="T")

        self.propeller = Propeller(aircraft="conventional",
                                   conditions=self.conditions, reference=self.reference, geometry=self.prop_geometry,
                                   power_condition=self.power_condition, propulsor_type="conventional",
                                   bem_input=self.bem_input)

        self.nacelle = Nacelle(geometry=self.nacelle_geometry, conditions=self.conditions, reference=self.reference,
                               n_blades=self.prop_geometry[0], prop_airfoil=self.prop_geometry[3])
    """ -------------------------------- geometric properties ------------------------------------------------------ """
    def area_wet(self):
        s_wet_htail = self.ht_tail.area_wetted()
        s_wet_vtail = self.vt_tail.area_wetted()

        s_wet_tot = s_wet_htail + s_wet_vtail
        return s_wet_tot

    """ ----------------------- vectors for plots ------------------------------------------------------------------ """
    def cd0_vector(self):
        cd0_nacelle = self.nacelle.cd0()[1]
        cd0_ht = self.ht_tail.cd0()[1]
        cd0_vt = self.vt_tail.cd0()[1]
        return [cd0_nacelle, cd0_ht, cd0_vt]

    def cd_interference_vector(self):
        cd_interference_ht = self.ht_tail.cd_interference()
        cd_interference_vt = self.vt_tail.cd_interference()
        cd_interference_nac = self.nacelle.cd_interference() * 2
        return [cd_interference_nac, cd_interference_ht, cd_interference_vt]

    """ ------------------------------- coefficient sums ------------------------------------------------------- """
    def cd_interference(self):
        cd_int_conv = self.ht_tail.cd_interference() + self.nacelle.cd_interference() + self.vt_tail.cd_interference()
        return cd_int_conv

    def cl_a(self):
        """ assume htail is only producing significant cl_a - > per radian"""
        cl_a_emp = self.ht_tail.cl_al()
        return cl_a_emp

    def cd_sum(self):
        """ regular cd-values """
        cd_htail = self.ht_tail.cd()[0]
        cd_vtail = self.vt_tail.cd()[0]
        return cd_htail + cd_vtail

    def ht_volume_coefficient(self):
        tv = tail_volume(self.ht_tail.area(), ref.lever_h, self.reference[0], ref.c_root_w)
        return tv

    def vt_volume_coefficient(self):
        tv = tail_volume(self.vt_tail.area(), ref.lever_h, self.reference[0], ref.b_w)
        return tv

    """ ------------------------------------------ Forces ------------------------------------------------------- """
    def thrust(self):
        thrust_conv = 2 * self.propeller.thrust()
        return thrust_conv

    def drag(self):
        drag_conv = self.ht_tail.drag_force() + self.vt_tail.drag_force()
        return drag_conv

    def lift(self):
        lift_conv = self.ht_tail.lift_force()
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

    def weight_fan(self):
        w_fan = 2 * self.propeller.weight_fan()
        return w_fan

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
    cl_norm = []
    cd_norm = []


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
                                    altitude=7000,
                                    bem_input=[41420.85603250924, 26482.06279555917, -1.4475057750305072e-13, 
                                               0.8892292886261024, 0.3297344029147765, 0.3201225439053968, 5, 10])

        cl.append(emp.ht_tail.cl())
        cd.append(emp.ht_tail.cd() + emp.vt_tail.cd())
        cl_norm.append(emp.ht_tail.cl_prime())
        cd_norm.append(emp.vt_tail.cd_prime() + emp.ht_tail.cd_prime())

    print(cl)
    print(cd)
    print(cl_norm)
    print(cd_norm)


    plt.figure('CL - alpha - Empennage comparison')
    plt.plot(a, cl, label=r'T-tail ATR72-600', color="tab:blue")
    plt.xlim([0, 20])
    # plt.ylim([0, 0.5])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(r'$C_{L}$ vs. $\alpha$')
    plt.legend()
    plt.grid(True)

    plt.figure('CD - alpha - Empennage comparison')
    plt.plot(a, cd, label=r'T-tail ATR72-600', color="tab:blue")
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.title(r'$C_{D}$ vs. $\alpha$')
    plt.legend()
    plt.grid(True)

    plt.figure('CL - CD - Empennage comparison')
    plt.plot(cd, cl, label=r'T-tail ATR72-600', color="tab:blue")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    # plt.xlim([0, 0.04])
    # plt.ylim([-0.1, 0.3])
    plt.title(r'$C_{L}$ vs. $C_{D}$')
    plt.legend()
    plt.grid(True)
    plt.show()
    """
