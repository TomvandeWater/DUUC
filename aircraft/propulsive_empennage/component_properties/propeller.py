import numpy as np
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa
from analysis_modules.factors import skin_friction
import data.atr_reference as ref
import matlab.engine
import matplotlib.pyplot as plt
# BEM_matlab_engine = matlab.engine.start_matlab()


class Propeller:
    """ This is the propeller class for the propulsive empennage"""
    def __init__(self, aircraft: str, conditions, reference, geometry, rpm: float, power_condition: str, va_inlet: float,
                 propulsor_type: str, bem_input, duct_diameter: float):
        super().__init__()
        self.aircraft = aircraft
        self.n_blades = geometry[0]
        self.prop_diameter = geometry[1]
        self.hub_diameter = geometry[2]
        self.prop_airfoil = geometry[3]
        self.prop_sweep = geometry[4]
        self.prop_pitch = geometry[5]
        self.c_root = geometry[6]
        self.c_tip = geometry[7]

        self.pc = power_condition
        self.propulsor_type = propulsor_type
        self.ref_area = reference[0]

        self.v_inf = conditions[0]
        self.alpha = conditions[1]
        self.altitude = conditions[2]
        self.mach = conditions[3]
        self.density = air_density_isa(self.altitude)[0]

        self.bem_input = bem_input
        self.T = self.bem_input[0]
        self.N = self.bem_input[6]
        self.duct_diameter = duct_diameter
        self.va_inlet = va_inlet
        self.rpm = rpm
    """ ---------------------------------- INFLOW PROPERTIES ------------------------------------------------------ """
    def inflow_velocity(self):
        if self.aircraft == "DUUC":
            if self.pc == "off":
                v_ind = (self.v_inf * np.pi * self.duct_diameter ** 2 * 0.25) / (np.pi / 4 * self.prop_diameter ** 2)
                u_prop = v_ind
                return u_prop
            else:
                v_ind = (self.v_inf * np.pi * self.duct_diameter ** 2 * 0.25) / (np.pi / 4 * self.prop_diameter ** 2)
                u_prop = v_ind
                return u_prop
        elif self.aircraft == "conventional":
            v_ind = self.v_inf
            return v_ind
        else:
            TypeError("Wrong Aircraft type specified")

    def inflow_angle(self):
        inflow_prop = self.alpha / 2

        angle_rad = np.radians(inflow_prop)
        return inflow_prop, angle_rad

    def reynolds_number(self):
        re_blade = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.c_root)
        return re_blade

    """ ----------------------------------- GEOMETRIC PROPERTIES --------------------------------------------------- """
    def area(self):
        area_prop = np.pi * (self.prop_diameter / 2) ** 2
        return area_prop

    def u_wake(self):
        den = 2 / self.density
        t_s = self.thrust() / self.area()
        u1_prop = np.sqrt(den * t_s + self.inflow_velocity() ** 2)
        return u1_prop

    @staticmethod
    def t_c():
        """ defined on comparable propeller blade"""
        t_c_prop_root = 0.10
        t_c_prop_tip = 0.08
        t_c_av = (t_c_prop_tip + t_c_prop_root) / 2
        return t_c_prop_root, t_c_prop_tip, t_c_av

    def area_wetted(self):
        wet_prop_bl = (2 * (1 + 0.5 * self.t_c()[2]) * self.prop_diameter / 2 * (self.c_root + self.c_tip) / 2)

        wet_prop = self.n_blades * wet_prop_bl
        return wet_prop

    """ ------------------------------ RATIOS USED FOR NORMALIZATION ----------------------------------------------- """
    def area_ratio(self):
        ar_prop = self.area() / self.ref_area
        return ar_prop

    def area_ratio_wet(self):
        ar_w_prop = self.area_wetted() / self.ref_area
        return ar_w_prop

    def velocity_ratio(self):
        v_ratio = self.inflow_velocity() ** 2 / self.v_inf ** 2
        return v_ratio

    """ ---------------------------------------- COEFFICIENTS ----------------------------------------------------- """
    def cn(self):
        cn_prop = self.N / (0.5 * air_density_isa(self.altitude)[0] * self.inflow_velocity() ** 2 * self.area())
        return cn_prop

    def ct(self):
        ct_prop = self.T / (0.5 * air_density_isa(self.altitude)[0] * self.inflow_velocity() ** 2 * self.area())
        return ct_prop

    def cd0(self):
        if self.pc == "off":
            """ in power off conditions the drag is estimated by only the skin friction drag"""
            cd0_prop = skin_friction(self.reynolds_number(), "t") * self.n_blades
            return cd0_prop
        else:
            """ in power on conditions the skin friction drag is assumed to be neglible"""
            cd0_prop = 0
            return cd0_prop

    def cd(self):
        inflow = self.inflow_angle()[1]

        if self.pc == "off":
            cd_prop = (self.cd0() * np.cos(inflow) * self.n_blades)
            cd_norm = cd_prop * self.area_ratio_wet() * self.velocity_ratio()
            return cd_prop, cd_norm
        else:
            """ installed drag (thrust)"""
            cd_prop = self.ct() * np.cos(inflow) - self.cn() * np.sin(inflow)
            cd_norm = cd_prop * self.area_ratio() * self.velocity_ratio()
            return cd_prop, cd_norm

    def cl(self):
        inflow = self.inflow_angle()[1]

        if self.pc == "off":
            cl_prop = self.cd()[0] * np.sin(inflow) * self.n_blades
            cl_norm = cl_prop * self.area_ratio_wet() * self.velocity_ratio()
            return cl_prop, cl_norm
        else:
            cl_1 = self.cn() * np.cos(inflow)

            cl_2 = self.ct() * np.sin(inflow)

            cl_prop = cl_1 + cl_2
            cl_norm = cl_prop * self.area_ratio_wet() * self.velocity_ratio()
            return cl_prop, cl_norm

    """ ---------------------------------- INTERFERENCE EFFECTS ---------------------------------------------------- """
    def cd_interference(self):
        if self.pc == "off":
            norm_area_nac = ((self.t_c()[0] * self.c_root) * self.c_root) / self.ref_area
            norm_area_tip = ((self.t_c()[1] * self.c_tip) * self.c_tip) / self.ref_area
            norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

            cd_int_nac = (self.n_blades * drag_interference(self.t_c()[0], "plane")
                          * norm_speed * norm_area_nac)

            cd_int_tip = (self.n_blades * drag_interference(self.t_c()[1], "plane")
                          * norm_speed * norm_area_tip)

            cd_int_prop = cd_int_nac + cd_int_tip
            return cd_int_prop
        else:
            """ Assume that interference drag is a neglible contribution when power is on"""
            cd_int_prop = 0
            return cd_int_prop

    """ ---------------------------------------- FORCES ------------------------------------------------------------ """
    def thrust(self):
        thrust_prop = self.T
        return thrust_prop

    def lift_force(self):
        l_prop = self.cl()[0] * (0.5 * air_density_isa(self.altitude)[0] * self.inflow_velocity() ** 2 * self.area())
        return l_prop

    def drag_force(self):
        d_prop = self.cd()[0] * (0.5 * air_density_isa(self.altitude)[0] * self.inflow_velocity() ** 2 * self.area())
        return d_prop

    """ ----------------------------------------- WEIGHT ----------------------------------------------------------- """
    def weight_engine(self):
        """ based on 1 engine"""
        me = ref.m_eng
        ke = 1.35  # for propeller driven aircraft with more than 1 engine
        kthr = 1.18  # accounts for reverse thrust
        if self.propulsor_type == "conventional":
            m_eng = ke * kthr * me
            return m_eng
        elif self.propulsor_type == "hybrid":

            m_bat = 300000 / 250  # 250 Wh/kg estimate
            m_gen = 800 / 1.5  # 800 kW generator -> 1.5 kg/kW estimate
            m_elec = m_bat / 3.2
            m_mot = 2 * (750 * 1.2)  # 750 kW per motor -> 1.2 kg/kW estimate
            m_cool = 600

            m_eng = m_bat + m_gen + m_elec + m_mot + m_cool
            return m_eng
        else:
            raise ValueError("PE propeller.py -> Propulsor type not correctly specified or included in the model")

    def weight_fan(self):
        m_fan = ref.m_blade * self.n_blades
        m_fan = m_fan * 1.10  # 10 percent mass added for spinner
        return m_fan


"""----- test section -----"""
"""
alfa = np.linspace(0, 20, 41)
cl = []
cd = []
cl_norm = []
cd_norm = []
cn = []
ct = []

for i in range(len(alfa)):
    duuc: Propeller = Propeller(n_blades=6, prop_diameter=3.6,
                                hub_diameter=0.2, prop_airfoil='ARAD8',
                                prop_sweep=0, prop_pitch=2, rpm=1200,
                                power_condition="on",
                                va_inlet=750,
                                alpha=alfa[i],
                                s_ref=62,
                                c_root=0.2,
                                c_tip=0.2,
                                v_inf=128,
                                normal=1900,
                                thrust=5000,
                                altitude=7000,
                                propulsor_type='conventional')

    cl.append(duuc.cl()[0])
    cl_norm.append(duuc.cl()[1])
    cd.append(duuc.cd()[0])
    cd_norm.append(duuc.cd()[1])
    ct.append(duuc.ct())
    cn.append(duuc.cn())

plt.figure('Propeller coefficients')
plt.plot(alfa, cl, label=r'Cl')
plt.plot(alfa, cd, label=r'Cd')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'Coefficients [-]')
plt.title(r'Propeller coefficients')
plt.legend()
plt.grid(True)

plt.figure('Propeller coefficients normalized')
plt.plot(alfa, cl_norm, label=r'Cl - norm')
plt.plot(alfa, cd_norm, label=r'Cd - norm')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'Coefficients [-]')
plt.title(r'Propeller coefficients')
plt.legend()
plt.grid(True)

plt.figure('Propeller coefficients 2')
plt.plot(alfa, cn, label=r'Cn')
plt.plot(alfa, ct, label=r'Ct')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'Coefficients [-]')
plt.title(r'Propeller coefficients 2')
plt.legend()
plt.grid(True)


plt.show() """
