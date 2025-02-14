import flow_conditions
from aircraft.propulsive_empennage.component_properties.pylon import Pylon
from aircraft.propulsive_empennage.component_properties.nacelle import Nacelle
from aircraft.propulsive_empennage.component_properties.support import SupportStrut
from aircraft.propulsive_empennage.component_properties.duct import Duct
from aircraft.propulsive_empennage.component_properties.control_vane import ControlVane
from aircraft.propulsive_empennage.component_properties.propeller import Propeller


class PropulsiveEmpennage:
    def __init__(propemp, rpm: float, alpha: float, power_condition: str, va_inlet: float,
                 re_duct: float, n_blades: float, prop_diameter: float, hub_diameter: float,
                 prop_airfoil: str, prop_sweep: float, prop_pitch: float, c_root: float,
                 c_tip: float, duct_diameter: float, duct_chord: float, duct_profile: str,
                 cant_angle: float, pylon_length: float, pylon_chord: float, pylon_profile: str,
                 d_exit: float, nacelle_length: float, nacelle_diameter: float,
                 support_length: float, support_chord: float, support_profile: str, hcv_span: float,
                 hcv_chord: float, control_profile: str, vcv_span: float, vcv_chord: float,
                 propulsor_type: str):
        super().__init__()
        # operating conditions
        propemp.rpm = rpm
        propemp.alpha = alpha
        propemp.power_condition = power_condition
        propemp.va_inlet = va_inlet
        propemp.re = re_duct
        propemp.cant = cant_angle
        propemp.d_exit = d_exit
        propemp.prop_type = propulsor_type

        # propeller properties
        propemp.n_bl = n_blades
        propemp.d_prop = prop_diameter
        propemp.d_hub = hub_diameter
        propemp.prop_airfoil = prop_airfoil
        propemp.s_prop = prop_sweep
        propemp.p_prop = prop_pitch
        propemp.c_root = c_root
        propemp.c_tip = c_tip

        # duct properties
        propemp.d_duct = duct_diameter
        propemp.c_duct = duct_chord
        propemp.airfoil_duct = duct_profile

        # pylon properties
        propemp.l_pylon = pylon_length
        propemp.c_pylon = pylon_chord
        propemp.pylon_profile = pylon_profile

        # nacelle properties
        propemp.l_nacelle = nacelle_length
        propemp.d_nacelle = nacelle_diameter

        # support properties
        propemp.l_support = support_length
        propemp.c_support = support_chord
        propemp.support_profile = support_profile

        # horizontal control vane properties
        propemp.b_hcv = hcv_span
        propemp.c_hcv = hcv_chord
        propemp.hcv_profile = control_profile

        # vertical control vane properties
        propemp.b_vcv = vcv_span
        propemp.c_vcv = vcv_chord
        propemp.vcv_profile = control_profile

        s_ref = propemp.d_duct * propemp.c_duct

        # Initiate propeller class
        propemp.propeller = Propeller(propemp.n_bl, propemp.d_prop, propemp.d_hub,
                                      propemp.prop_airfoil, propemp.s_prop, propemp.p_prop,
                                      propemp.rpm, propemp.power_condition, propemp.va_inlet,
                                      propemp.alpha, propemp.re, s_ref, propemp.c_root,
                                      propemp.c_tip)

        # calculate properties based on input from propeller
        u1 = propemp.propeller.u1()
        tc_prop = propemp.propeller.tc()
        v_prop = propemp.propeller.inflow_velocity()
        cn_prop = propemp.propeller.cn()
        u_mom = 0.5 * (v_prop + u1)

        # Initiate duct class
        propemp.duct = Duct(propemp.d_duct, propemp.c_duct, propemp.airfoil_duct, propemp.alpha,
                            propemp.re, propemp.power_condition, u_mom, tc_prop)

        # calculate properties based on input from duct
        ref_area = propemp.duct.proj_area()

        # Initiate pylon class
        propemp.pylon = Pylon(propemp.l_pylon, propemp.c_pylon, propemp.pylon_profile,
                              propemp.power_condition, propemp.cant, propemp.alpha, ref_area)

        # Initiate nacelle class
        propemp.nacelle = Nacelle(propemp.l_nacelle, propemp.d_nacelle, propemp.prop_type,
                                  propemp.re, propemp.power_condition, u_mom, alpha, ref_area)

        # Initiate support class
        propemp.support = SupportStrut(propemp.l_support, propemp.c_support,
                                       propemp.support_profile, propemp.cant,
                                       propemp.power_condition, v_prop, u_mom, alpha, tc_prop,
                                       cn_prop, ref_area)

        # Initiate control vane class for elevator (1 piece)
        propemp.elevator = ControlVane(propemp.b_hcv, propemp.c_hcv, propemp.hcv_profile,
                                       propemp.power_condition, propemp.va_inlet,
                                       propemp.d_exit, u1, ref_area, flow_conditions.delta_e,
                                       propemp.re)

        # Initiate control vane class for rudder (1 piece)
        propemp.rudder = ControlVane(propemp.b_vcv, propemp.c_vcv, propemp.vcv_profile,
                                     propemp.power_condition, propemp.va_inlet, propemp.d_exit, u1,
                                     ref_area, flow_conditions.delta_r, propemp.re)

    def cl_prime(self):
        cl_prime_pe = (self.pylon.cl_prime() + self.duct.cl_prime() + self.support.cl_prime()
                       + self.nacelle.cl_prime() + 2 * self.elevator.cl_prime())
        return cl_prime_pe

    def cd_prime(self):
        cd_prime_pe = (self.pylon.cd_prime() + self.duct.cd_prime() + self.support.cd_prime()
                       + self.nacelle.cd_prime() + 2 * self.elevator.cd_prime()
                       + 2 * self.rudder.cd_prime() + self.propeller.cd())

        cd_interference = (self.pylon.cd_interference() + self.support.cd_interference()
                           + self.propeller.cd_interference() + 2 * self.rudder.cd_interference()
                           + 2 * self.elevator.cd_interference())

        cd_total = cd_prime_pe + cd_interference

        return cd_total
