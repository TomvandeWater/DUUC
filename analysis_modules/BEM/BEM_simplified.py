import numpy as np
import os
import data.atr_reference as ref
import matplotlib.pyplot as plt


class BEM:
    def __init__(self, radius, num_blades, rpm, prop_airfoil):
        self.R = radius
        self.B = num_blades
        self.omega = rpm * 2 * np.pi / 60
        self.prop_airfoil = prop_airfoil

        # Simplified blade geometry (can be customized)
        self.chord = lambda r: 0.1 * (1 - r / self.R) + 0.125
        self.twist = lambda r: 20 * (1 - r / self.R)

    def calculate_thrust(self, V_inf, rho=0.589, num_elements=50):
        r = np.linspace(0.1 * self.R, self.R, num_elements)
        dr = r[1] - r[0]

        thrust = 0
        torque = 0
        normal_f = 0

        for i in range(num_elements):
            # Local blade geometry
            chord = self.chord(r[i])
            twist = np.radians(self.twist(r[i]))

            # Initial guess for induction factors
            a = 0.1
            a_prime = 0.01

            # Iterate to find induction factors
            for _ in range(100):
                phi = np.arctan2(V_inf * (1 + a), self.omega * r[i] * (1 - a_prime))
                alpha = phi - twist
                # print(f"Alpha: {alpha}")
                Cl, Cd = self.get_airfoil_data(alpha)
                # print(f"Cl: {Cl}, Cd: {Cd}")
                Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
                Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
                # print(f"Cn: {Cn}")
                sigma = self.B * chord / (2 * np.pi * r[i])
                F = self.tip_loss_correction(r[i], phi)

                a_new = 1 / ((4 * F * np.sin(phi) ** 2) / (sigma * Cn) + 1)
                a_prime_new = 1 / ((4 * F * np.sin(phi) * np.cos(phi)) / (sigma * Ct) - 1)

                if abs(a - a_new) < 1e-5 and abs(a_prime - a_prime_new) < 1e-5:
                    break

                a, a_prime = a_new, a_prime_new

            # Calculate local thrust and torque
            W = np.sqrt((V_inf * (1 + a)) ** 2 + (self.omega * r[i] * (1 - a_prime)) ** 2)
            dT = 0.5 * rho * W ** 2 * self.B * chord * Cn * dr
            dQ = 0.5 * rho * W ** 2 * self.B * chord * Ct * r[i] * dr

            thrust += dT
            torque += dQ
            normal_f = torque / self.R

        return thrust, torque, normal_f

    def tip_loss_correction(self, r, phi):
        f = self.B / 2 * (self.R - r) / (r * np.sin(phi))
        # Calculate cos_value
        cos_value = np.exp(-f)

        # Clamp cos_value to the range [-1, 1]
        cos_value = np.clip(cos_value, -1, 1)

        return 2 / np.pi * np.arccos(cos_value)

    def get_airfoil_data(self, alpha):
        filepath = os.path.join(r"/data/Polars",
                                f"{self.prop_airfoil}_polar.txt")
        # Normalize the filepath to ensure it resolves correctly
        filepath = os.path.abspath(filepath)

        data = np.loadtxt(filepath, skiprows=2)
        polar_alpha = data[:, 0]
        polar_cl = data[:, 1]
        polar_cd = data[:, 2]

        # Interpolate airfoil data
        cl = np.interp(np.degrees(alpha), polar_alpha, polar_cl)
        cd = np.interp(np.degrees(alpha), polar_cd, polar_cd)
        return cl, cd

"""
v_cruise = 128  # cruise speed [m/s] circa 250 knts
v_array = np.linspace(v_cruise - 20, v_cruise + 20, 41)
thrust_array = []

for i in range(len(v_array)):
    u_inf = float(v_array[i])
    prop = BEM(radius=ref.blade_diameter/2, num_blades=ref.n_blades, rpm=1000, prop_airfoil="Hamilton568F")

    thrust, torque, normal_f = prop.calculate_thrust(u_inf)

    thrust_array.append([u_inf, thrust])

plot_array = np.array(thrust_array)

plt.figure('Thrust plot')
plt.plot(plot_array[:, 0], plot_array[:, 1], label=r'Total drag coefficient')
plt.xlabel(r'$v_{inf}$')
plt.ylabel(r'$Thrust [N]$')
plt.title(r'Thrust array')
plt.legend()
plt.grid(True)
plt.show()"""


# Example usage
if __name__ == "__main__":
    prop = BEM(radius=1.8, num_blades=6, rpm=1000, prop_airfoil="Hamilton568F")

    V_inf = 128  # m/s
    thrust, torque, normal_f = prop.calculate_thrust(V_inf)

    print(f"Thrust: {thrust:.2f} N")
    print(f"Torque: {torque:.2f} N·m")
    print(f"Power: {torque * prop.omega:.2f} W")
    print(f"Normal: {normal_f:.2f}")


