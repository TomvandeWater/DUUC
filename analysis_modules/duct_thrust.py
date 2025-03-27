import numpy as np
import matplotlib.pyplot as plt

# Constants
rho = 1.225  # Air density at sea level (kg/m^3)
D_p = 3.5  # Propeller diameter (m)
n = 1000  # Rotational speed (rpm)
V_0 = 80  # Cruise speed (m/s)


# Function to calculate thrust components
def calculate_thrust(D_duct, L_duct, V_0, J, delta):
    # Propeller thrust
    A_p = np.pi * (D_p / 2) ** 2
    V_p = n * np.pi * D_p / 60  # Propeller tip speed
    V_e = V_p + V_0  # Simplified exit velocity

    # Limit V_e to a realistic value to avoid excessive thrust
    V_e = min(V_e, 200)  # m/s

    T_p = rho * A_p * V_p * (V_e - V_0)

    # Duct pressure thrust
    A_exit = np.pi * (delta * D_duct / 2) ** 2
    T_duct_pressure = 0.5 * rho * A_exit * (V_e ** 2 - V_0 ** 2)

    # Duct shape thrust (revised to be more realistic)
    T_duct_shape = 0.001 * rho * n ** 2 * D_duct ** 4

    # Total thrust
    T_total = T_p + T_duct_pressure + T_duct_shape

    return T_total, T_duct_pressure, T_duct_shape


# Sensitivity analysis
D_ducts = np.linspace(2.5, 4.5, 100)  # Duct diameters
deltas = np.linspace(0.7, 1.2, 100)  # Expansion ratios

# Duct diameter sensitivity
C_T_duct_diameter = []
C_T_pressure_contribution = []
C_T_shape_contribution = []

for D_duct in D_ducts:
    T_total, T_pressure, T_shape = calculate_thrust(D_duct, 2.0, V_0, 0.3, 0.9)

    # Dynamic pressure
    q = 0.5 * rho * V_0 ** 2

    # Reference area for coefficient calculation
    A_ref = np.pi * (D_duct / 2) ** 2  # Use duct projected area as reference

    # Thrust coefficients
    C_T_total = T_total / (q * A_ref)
    C_T_pressure = T_pressure / (q * A_ref)
    C_T_shape = T_shape / (q * A_ref)

    C_T_duct_diameter.append(C_T_total)
    C_T_pressure_contribution.append(C_T_pressure)
    C_T_shape_contribution.append(C_T_shape)

# Expansion ratio sensitivity
C_T_expansion_ratio = []
C_T_expansion_pressure_contribution = []
C_T_expansion_shape_contribution = []

for delta in deltas:
    T_total, T_pressure, T_shape = calculate_thrust(3.0, 2.0, V_0, 0.3, delta)

    # Dynamic pressure
    q = 0.5 * rho * V_0 ** 2

    # Reference area for coefficient calculation
    A_ref = np.pi * (3.0 / 2) ** 2  # Use duct projected area as reference

    # Thrust coefficients
    C_T_total = T_total / (q * A_ref)
    C_T_pressure = T_pressure / (q * A_ref)
    C_T_shape = T_shape / (q * A_ref)

    C_T_expansion_ratio.append(C_T_total)
    C_T_expansion_pressure_contribution.append(C_T_pressure)
    C_T_expansion_shape_contribution.append(C_T_shape)

# Plot figures
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(D_ducts, C_T_duct_diameter, label='Total Thrust Coefficient')
plt.plot(D_ducts, C_T_pressure_contribution, linestyle='--', label='Pressure Contribution Coefficient')
plt.plot(D_ducts, C_T_shape_contribution, linestyle='--', label='Shape Contribution Coefficient')
plt.xlabel('Duct Diameter (m)')
plt.ylabel('Thrust Coefficient')
plt.title('Thrust Sensitivity to Duct Diameter')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(deltas, C_T_expansion_ratio, label='Total Thrust Coefficient')
plt.plot(deltas, C_T_expansion_pressure_contribution, linestyle='--', label='Pressure Contribution Coefficient')
plt.plot(deltas, C_T_expansion_shape_contribution, linestyle='--', label='Shape Contribution Coefficient')
plt.xlabel('Expansion Ratio')
plt.ylabel('Thrust Coefficient')
plt.title('Thrust Sensitivity to Expansion Ratio')
plt.legend()

plt.tight_layout()

plt.show()


def cT_duct_pressure(v_in, rpm, diameter, rho):
    n = rpm / 60
    cT_d_list = []
    for v in v_in:
        j = v / (n * diameter)
        if j < 0.35:
            T_d = 0.26 * rho * rpm ** 2 * diameter ** 4
            cT_d = T_d / (0.5 * rho * v ** 2 * 62)
        else:
            T_d = -0.04 * rho * rpm ** 2 * diameter ** 4
            cT_d = T_d / (0.5 * rho * v ** 2 * 62)
        cT_d_list.append(cT_d)
    return np.array(cT_d_list)


v_inf = np.linspace(0.1, 120.1, 121)

plt.figure(figsize=(8, 6))
plt.plot(v_inf, cT_duct_pressure(v_inf, 1000, 3.8, 1.225), label='RPM = 1000')
plt.plot(v_inf, cT_duct_pressure(v_inf, 2000, 3.8, 1.225), label='RPM = 2000')
plt.plot(v_inf, cT_duct_pressure(v_inf, 3000, 3.8, 1.225), label='RPM = 3000')
plt.plot(v_inf, cT_duct_pressure(v_inf, 4000, 3.8, 1.225), label='RPM = 4000')
plt.xlabel(r'$V_{\infty}$ [m/s]')
plt.ylabel(r'$C_{Td}$ [-]')
plt.title('Duct thrust coefficient due to pressure')
plt.legend()
plt.show()


