from analysis_modules.aerodynamic import c_tip, c_root, drag_zero_lift
""" Reference data from the fact sheet ATR72-600 and from this 
website: https://www.fzt.haw-hamburg.de/pers/Scholz/arbeiten/DataHirsch/Aircraft_Database.html"""
import numpy as np

""" Total dimensions"""
l_ac = 27.17  # length aircraft? [m]
b_w = 27.05  # wing span [m]
h_ac = 7.65  # height aircraft [m]
l_cab = 19.21  # cabin length [m]
diameter = 2.57  # diameter of the fuselage [m]

""" Wing parameters"""
s_w = 61  # wing area [m^2]
tr_w = 0.47  # taper ratio wing [-]
phi_qc_w = 3.5  # quarter chord sweep [deg]
c_root_w = c_root(s_w / 2, b_w / 2, tr_w)
c_tip_w = c_tip(c_root_w, tr_w)
ar_w = b_w ** 2 / s_w
cl_wing = -0.116  # clean wing lift coefficient M = 0.0792 Re 2.005e5
alpha_install_wing = 4  # installation angle of the wing (w.r.t zero lift line) onto the fuselage [deg]
wing_airfoil = "43013"  # wing airfoil profile from the Naca 5 series -> or 43018

""" Horizontal tail surface"""
s_ht = 11.7  # horizontal tail surface [m^2]
# s_ht = 10.7  # from Hamburg university research
b_h = 7.31  # span horizontal tail surface [m]
tr_h = 0.525  # taper ratio [-]
phi_qc_h = 5.5  # quarter chord sweep [deg]
phi_hc_h = 25  # half chord sweep [deg]
ar_h = 6  # aspect ratio horizontal stab [-]
c_root_h = c_root(s_ht, b_h, tr_h)
c_tip_h = c_tip(c_root_h, tr_h)
tail_volume_h = 1.05  # tail volume coefficient horizontal tail (from Hamburg University)
airfoil_ht = "0009"  # horizontal tailplane airfoil profile from NACA 4 series
lever_h = 13.565  # leverage arm between tail [m]

""" Vertical tail surface"""
s_vt = 14.9  # vertical tail surface [m^2]
# s_vt = 14.90  # from Hamburg university research
b_v = 4.43  # tail height [m]
tr_v = 0.65  # taper ratio [-]
phi_qc_v = 29  # quarter chord sweep [deg]
phi_hc_v = 24  # half chord sweep [deg]
ar_v = 1.6  # aspect ratio vertical stab [-]
c_root_v = c_root(s_vt, b_v, tr_v)
c_tip_v = c_tip(c_root_v, tr_v)
tail_volume_v = 0.119  # tail volume coefficient vertical tail (from Hamburg University)
airfoil_vt = "0012"  # vertical tailplane airfoil profile from NACA 4 series

""" Fuselage"""
h_f = 2.63  # height of the fuselage [m]
w_f = 2.87  # width of the fuselage [m]
diameter_fuselage = 2.77  # fuselage diameter [m]
l_cabin = 19.25  # cabin length [m]
l_cockpit = 3.880  # cockpit length [m] from raydome to front door
l_tail = 8.310  # tail length [m]

""" Engines PW Canada PW127 M/N"""
P_TO = 2475  # Take off power [SHP]
P_TO_EO = 2750  # Take off power one engine [SHP]
P_max_con = 2500  # Power max continuous [SHP]
P_max_climb = 2192  # Power max climb [SHP]
P_max_cruise = 2132  # Power max cruise [SHP]
l_nacelle = 3.4425  # nacelle length [m] (based on dimensional drawing)
d_nacelle = 1.1475  # nacelle diameter (average) [m] (based on dimensional drawing)
                    # -> assume constant diameter
m_eng = 481  # engine mass in [kg]

""" Blades Hamilton Standard 568F """
blade_diameter = 3.93  # blade diameter [m]
n_blades = 6  # number of blades [-]
propeller_airfoil = "Hamilton568F"  # airfoil propeller blade
m_blade = 52 * 0.45359237  # blade mass [kg]

""" Weights """
MTOW = 22800  # Maximum take-off weight (basic) [kg]
ZFW_max = 20800  # Max zero fuel weight (basic) [kg]
OEW = 13010  # Operational Empty weight (Tech spec.) [kg]
payload_max = 8550  # Max payload (at typical in-service OEW) [kg]
load_fuel = 5000  # Max fuel load [kg]

""" En-route performance"""
v_cruise = 510  # Max cruise speed (95% MTOW - ISA - Optimum FL) [km/h]
m_fuel = 762  # Fuel flow at cruise speed [kg/hr]
range_reg = 758  # Range with max pax [NM]
v_dive = 130  # dive speed [m/s]

""""  Weight estimations based on Hamburg university """
m_wing = 3045  # wing mass [kg]
m_fuse = 2323  # fuselage mass [kg]
m_ht = 124  # horizontal tail mass [kg]
m_vt = 179  # vertical tail mass [kg]
m_lg = 961  # landing gear mass [kg]
m_nc = 242  # engine nacelle mass [kg]
m_en = 1533  # installed engine mass [kg]
m_sys = 3114  # systems mass [kg]
m_sup = 1050  # supplemental mass [kg]

""" Drag estimations"""
drag_fuse = [2.24e-3, 1.088, 1, 3.3]
drag_wing = [3.56e-3, 1.84, 1, 2.08]
drag_ht = [3.392e-3, 1.368, 1.04, 0.17]
drag_vt = [3.933e-3, 1.419, 1.04, 0.22]
drag_nac = [3.292e-3, 1.072, 1.5, 0.3]

fuse_cd0 = drag_zero_lift(drag_fuse[0], drag_fuse[1], drag_fuse[2], drag_fuse[3])
wing_cd0 = drag_zero_lift(drag_wing[0], drag_wing[1], drag_wing[2], drag_wing[3])
ht_cd0 = drag_zero_lift(drag_ht[0], drag_ht[1], drag_ht[2], drag_ht[3])
vt_cd0 = drag_zero_lift(drag_vt[0], drag_vt[1], drag_vt[2], drag_vt[3])
nac_cd0 = drag_zero_lift(drag_nac[0], drag_nac[1], drag_nac[2], drag_nac[3])


"""" Create reference matrix for AVL purposes"""
matrix = [airfoil_ht, airfoil_vt, c_root_h, c_root_v, tr_h, tr_v, phi_qc_h, phi_qc_v, b_h, b_v,
          m_ht, m_vt, ht_cd0+vt_cd0]
