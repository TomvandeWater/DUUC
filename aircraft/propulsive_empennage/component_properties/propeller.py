import numpy as np
from analysis_modules.aerodynamic import drag_interference, reynolds
from analysis_modules.ISA import air_density_isa
from analysis_modules.factors import skin_friction
import data.atr_reference as ref
import matlab.engine
BEM_matlab_engine = matlab.engine.start_matlab()


class Propeller:
    """ This is the propeller class for the propulsive empennage"""
    def __init__(self, n_blades: float, prop_diameter: float, hub_diameter: float,
                 prop_airfoil: str, prop_sweep: float, prop_pitch: float, rpm: float,
                 power_condition: str, va_inlet: float, alpha: float,
                 s_ref: float, c_root: float, c_tip: float, v_inf: float, propulsor_type: str, altitude: float):
        super().__init__()
        self.n_blades = n_blades
        self.prop_diameter = prop_diameter
        self.hub_diameter = hub_diameter
        self.prop_airfoil = prop_airfoil
        self.prop_sweep = prop_sweep
        self.prop_pitch = prop_pitch
        self.rpm = rpm
        self.pc = power_condition
        self.va_inlet = va_inlet
        self.alpha = alpha
        self.ref_area = s_ref
        self.c_root = c_root
        self.c_tip = c_tip
        self.v_inf = v_inf
        self.propulsor_type = propulsor_type
        self.altitude = altitude
        self.density = air_density_isa(self.altitude)[0]

    def inflow_velocity(self):
        if self.pc == "off":
            v_ind = self.va_inlet / (np.pi / 4 * self.prop_diameter ** 2)
            u_prop = v_ind
            return u_prop
        else:
            v_ind = self.va_inlet / (np.pi / 4 * self.prop_diameter ** 2)
            u_prop = v_ind
            return u_prop

    def inflow_angle(self):
        inflow_prop = self.alpha / 2
        return inflow_prop

    def reynolds_number(self):
        re_blade = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.c_root)
        return re_blade

    def area(self):
        area_prop = np.pi * (self.prop_diameter / 2) ** 2
        return area_prop

    @staticmethod
    def t_c():
        """ defined on comparable propeller blade"""
        t_c_prop_root = 0.10
        t_c_prop_tip = 0.08
        t_c_av = (t_c_prop_tip + t_c_prop_root) / 2
        return t_c_prop_root, t_c_prop_tip, t_c_av

    def wet_area(self):
        wet_prop_bl = (2 * (1 + 0.5 * self.t_c()[2]) * self.prop_diameter / 2
                       * (self.c_root + self.c_tip) / 2)

        wet_prop = self.n_blades * wet_prop_bl
        return wet_prop

    """ ----------------------- run BEM model -----------------------------------------------------------------------"""
    def BEM(self):
        #BEM_matlab_engine.cd(r'C:\Users\tomva\pythonProject\DUUC\analysis_modules\BEM')
        #T_out, Q_out, N_out, Tc, Cp, CT = BEM_matlab_engine.BEM2(6, 3.9, 10, 0, 50, 0, 0.6, 1.225, 15, 'HM568F',
                                                                 #nargout=6)

        #BEM_matlab_engine.quit()

        T_out =1
        Q_out =1
        N_out = 1
        Tc = 1
        Cp = 1
        CT = 1
        return T_out, Q_out, N_out, Tc, Cp, CT

    """ Coefficient calculations """
    def cd0(self):
        if self.pc == "off":
            """ in power off conditions the drag is estimated by only the skin friction drag"""
            cd0_prop = skin_friction(self.reynolds_number(), "t") * self.n_blades
            return cd0_prop
        else:
            """ in power on conditions the skin friction drag is assumed to be neglible"""
            cd0_prop = 0
            return cd0_prop

    def cd_prime(self):
        a = np.radians(self.alpha)
        inflow_a = np.radians(self.inflow_angle())
        norm_area = self.wet_area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        if self.pc == "off":
            cd_prime = (self.cd0() * np.cos(a - inflow_a) * self.n_blades * norm_area
                        * norm_speed)
            return cd_prime
        else:

            cd_prime = self.cn() * np.sin(a - inflow_a) * norm_area * norm_speed
            return cd_prime

    def cl_prime(self):
        a = np.radians(self.alpha)
        inflow_a = np.radians(self.inflow_angle())
        norm_area = self.wet_area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        if self.pc == "off":
            cl_prime = (self.cd0() * np.sin(a - inflow_a) * self.n_blades * norm_area
                        * norm_speed)
            return cl_prime
        else:
            cl_1 = self.tc() * np.sin(a)

            cl_2 = self.cn() * np.cos(a - inflow_a) * norm_area * norm_speed

            cl_prime = cl_1 + cl_2
            return cl_prime

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

    def thrust(self):
        """
        prop = BEM(radius=self.prop_diameter/2,
                   num_blades=self.n_blades,
                   rpm=self.rpm,
                   prop_airfoil=self.prop_airfoil)
        thrust_prop, torque_prop, normal_f = prop.calculate_thrust(self.v_inf)"""
        thrust_prop = 8000
        return thrust_prop

    def tc(self):
        tc_prop = self.thrust() / (0.5 * self.density * self.inflow_velocity() ** 2 * self.area())
        return tc_prop

    def cn(self):
        """
        prop = BEM(radius=self.prop_diameter / 2,
                   num_blades=self.n_blades,
                   rpm=self.rpm,
                   prop_airfoil=self.prop_airfoil)
        thrust_prop, torque_prop, normal_f = prop.calculate_thrust(self.v_inf)"""

        normal_f = 3200
        cn_prop = normal_f / (0.5 * self.density * self.inflow_velocity() ** 2
                              * self.area())
        return cn_prop

    def u1(self):
        u1_prop = np.sqrt((2/self.density) * (self.thrust() / self.area()) +
                          self.inflow_velocity() ** 2)
        return u1_prop

    """ Weight definition"""
    def weight_engine(self):
        """ based on 1 engine"""
        me = ref.m_eng
        ke = 1.35  # for propeller driven aircraft with more than 1 engine
        kthr = 1.18  # accounts for reverse thrust
        if self.propulsor_type == "conventional":
            m_eng = ke * kthr * me
            return m_eng
        if self.propulsor_type == "hybrid":
            m_eng = 10000000

            return m_eng
        else:
            print("PE propeller.py -> Propulsor type not correctly specified or included in the model")
            return None

    def weight_fan(self):
        m_fan = ref.m_blade * self.n_blades
        m_fan = m_fan * 1.10  # 10 percent mass added for spinner
        return m_fan


"""----- test section -----"""
"""
duuc: Propeller = Propeller(n_blades=3, prop_diameter=3.6,
                            hub_diameter=0.2, prop_airfoil='ARAD8',
                            prop_sweep=0, prop_pitch=2, rpm=6000,
                            power_condition="on",
                            va_inlet=750,
                            alpha=0,
                            s_ref=62,
                            c_root=0.2,
                            c_tip=0.2,
                            v_inf=30)

print(f"BEM output: {duuc.BEM()}") """
