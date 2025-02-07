import flow_conditions
import config
from analysis_modules.aerodynamic import reynolds
from analysis_modules.ISA import air_density_isa
from analysis_modules.xfoil_run import xfoil_polar
from propulsive_empennage.propulsive_empennage import PropulsiveEmpennage

""" Initialization of the analysis module"""
alpha1 = -5
alpha2 = 10
delta_u = 100

alpha = 0

""" Set operating conditions for several components"""
pylon_reynolds = reynolds(air_density_isa(flow_conditions.altitude),
                          flow_conditions.u_inf, config.pylon_chord, "Pylon")
xfoil_polar(config.pylon_airfoil, alpha1, alpha2, 0.25, pylon_reynolds,
            "pylon0012", flow_conditions.Mach)


""" setup the DUUC object"""
duuc: PropulsiveEmpennage = PropulsiveEmpennage(pylon_length=config.pylon_length,
                                                pylon_chord=config.pylon_chord,
                                                pylon_profile="Naca0006",
                                                cant_angle=config.cant_angle,
                                                alpha=alpha)
print(f"\n----- Pylon Forces -----")
duuc.pylon.lift()
duuc.pylon.drag()
duuc.pylon.moment()
duuc.pylon.thrust()
