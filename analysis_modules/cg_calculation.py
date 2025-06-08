import numpy as np
import matplotlib.pyplot as plt
from analysis_modules.plotting_functions import print_cg_mass
import data.atr_reference as ref
import os
import matplotlib.gridspec as gridspec
import config
import data.atr_reference as ref


class CenterOfGravity:
    """ determine center of gravity based on method of Scholz 2008, note: input of weight in [kg] and distance in [m]"""
    def __init__(self, w_lg_nose: float, w_lg_main: float, w_eng_wing: float, w_wing: float, w_fuse: float,
                 w_nac_wing: float, w_vt: float, w_ht: float, w_duct: float, w_sys: float, l_fuselage: float,
                 aircraft_type: str, c_mac_wing: float, x_duct: float, x_wing: float, w_fuel: float, w_pax: float,
                 z_PE: float):
        self.w_lg_nose = w_lg_nose
        self.w_lg_main = w_lg_main
        self.w_eng_wing = w_eng_wing
        self.w_wing = w_wing
        self.w_fuse = w_fuse
        self.w_nac_wing = w_nac_wing
        self.w_vt = w_vt
        self.w_ht = w_ht
        self.w_duct = w_duct
        self.aircraft_type = aircraft_type
        self.w_sys = w_sys
        self.l_fuselage = l_fuselage
        self.c_mac_wing = c_mac_wing
        self.x_duct = x_duct
        self.x_wing = x_wing
        self.w_fuel = w_fuel
        self.w_pax = w_pax
        self.z_PE = z_PE

    """ ---------------------------- determine center of gravity of each component -------------------------------- """
    def cg_loc_wing(self):
        if self.aircraft_type == "conventional":
            """ Data extracted from reference values ATR72-600"""
            cg_wing = 12.24
            cg_main = 12.429
            cg_eng_wing = 10.12
            cg_nac_wing = 10.12
            cg_fuel = 12.24
            return [cg_wing, cg_main, cg_eng_wing, cg_nac_wing, cg_fuel]
        if self.aircraft_type == "DUUC":
            cg_wing = self.x_wing + (self.c_mac_wing / 2)
            cg_main = cg_wing + 0.25
            cg_fuel = cg_wing
            cg_trim = 5
            return [cg_wing, cg_main, cg_fuel, cg_trim]
        else:
            return None

    def cg_loc_fus(self):
        if self.aircraft_type == "conventional":
            cg_fus = 0.45 * self.l_fuselage
            cg_vt = 23.99
            cg_ht = 23.99
            cg_nose = 1.75
            cg_sys = 11.5
            cg_pax = cg_fus
            return [cg_fus, cg_vt, cg_ht, cg_nose, cg_sys, cg_pax]
        if self.aircraft_type == "DUUC":
            cg_fus = 0.45 * self.l_fuselage
            cg_nose = 1.75
            cg_sys = 0.55 * self.l_fuselage
            cg_duct = self.x_duct
            cg_pax = cg_fus
            return [cg_fus, cg_nose, cg_sys, cg_duct, cg_pax]

    def cg_wing_group(self):
        if self.aircraft_type == 'conventional':
            mass_vect = [self.w_wing, self.w_lg_main, self.w_eng_wing * 2, self.w_nac_wing * 2, self.w_fuel]

            sum_w = sum(mass_vect)

            mcg_wing = self.w_wing * self.cg_loc_wing()[0]
            mcg_main = self.w_lg_main * self.cg_loc_wing()[1]
            mcg_eng_wing = self.w_eng_wing * 2 * self.cg_loc_wing()[2]
            mcg_nac_wing = self.w_nac_wing * 2 * self.cg_loc_wing()[3]
            mcg_fuel = self.w_fuel * self.cg_loc_wing()[4]

            mcg_vect = [mcg_wing, mcg_main, mcg_eng_wing, mcg_nac_wing, mcg_fuel]

            sum_mcg = sum(mcg_vect)

            cg_wing_group = sum_mcg / sum_w

            components = ["Wing", "Main", "Engine", "Nacelle", "Fuel"]
            file_path = r"C:\Users\tomva\pythonProject\DUUC\data\CG_breakdown"
            os.makedirs(file_path, exist_ok=True)
            input_file_path = os.path.join(file_path, f"ATR_wing_group_cg_breakdown.txt")
            if os.path.exists(input_file_path):
                os.remove(input_file_path)
            with open(input_file_path, "w") as file:
                file.write(print_cg_mass(components, mass_vect, self.cg_loc_wing(), mcg_vect, sum_w, sum_mcg, cg_wing_group,
                           "wing    ", "Conventional"))
                file.close()
            return cg_wing_group, sum_w

        if self.aircraft_type == "DUUC":
            mass_vect = [self.w_wing, self.w_lg_main, self.w_fuel]

            sum_w = sum(mass_vect)

            mcg_wing = self.w_wing * self.cg_loc_wing()[0]
            mcg_main = self.w_lg_main * self.cg_loc_wing()[1]
            mcg_fuel = self.w_fuel * self.cg_loc_wing()[2]

            mcg_vect = [mcg_wing, mcg_main, mcg_fuel]

            sum_mcg = sum(mcg_vect)

            cg_wing_group = sum_mcg / sum_w

            components = ["Wing", "Main", "Fuel"]
            file_path = r"C:\Users\tomva\pythonProject\DUUC\data\CG_breakdown"
            os.makedirs(file_path, exist_ok=True)
            input_file_path = os.path.join(file_path, f"DUUC_wing_group_cg_breakdown.txt")
            if os.path.exists(input_file_path):
                os.remove(input_file_path)
            with open(input_file_path, "w") as file:
                file.write(print_cg_mass(components, mass_vect, self.cg_loc_wing(), mcg_vect, sum_w, sum_mcg, cg_wing_group,
                           "wing     ", "DUUC   "))
                file.close()
            return cg_wing_group, sum_w
        else:
            return None

    def cg_fuselage_group(self):
        if self.aircraft_type == 'conventional':
            mass_vect = [self.w_fuse, self.w_vt, self.w_ht, self.w_lg_nose, self.w_sys, self.w_pax]

            sum_w = sum(mass_vect)

            mcg_fus = self.w_fuse * self.cg_loc_fus()[0]
            mcg_vt = self.w_vt * self.cg_loc_fus()[1]
            mcg_ht = self.w_ht * self.cg_loc_fus()[2]
            mcg_nose = self.w_lg_nose * self.cg_loc_fus()[3]
            mcg_sys = self.w_sys * self.cg_loc_fus()[4]
            mcg_pax = self.w_pax * self.cg_loc_fus()[5]

            mcg_vect = [mcg_fus, mcg_vt, mcg_ht, mcg_nose, mcg_sys, mcg_pax]

            sum_mcg = sum(mcg_vect)

            cg_fuse_group = sum_mcg / sum_w

            components = ["Fuselage", "Vertical Tail", "Horizontal Tail", "Nose", "Systems", "Pax"]

            file_path = r"C:\Users\tomva\pythonProject\DUUC\data\CG_breakdown"
            os.makedirs(file_path, exist_ok=True)
            input_file_path = os.path.join(file_path, f"ATR_fuselage_group_cg_breakdown.txt")
            if os.path.exists(input_file_path):
                os.remove(input_file_path)
            with open(input_file_path, "w") as file:
                file.write(print_cg_mass(components, mass_vect, self.cg_loc_fus(), mcg_vect, sum_w, sum_mcg, cg_fuse_group,
                           "fuselage", "Conventional"))
                file.close()
            return cg_fuse_group, sum_w

        if self.aircraft_type == "DUUC":
            mass_vect = [self.w_fuse, self.w_lg_nose, self.w_sys, self.w_duct, self.w_pax]
            sum_w = sum(mass_vect)

            mcg_fus = self.w_fuse * self.cg_loc_fus()[0]
            mcg_nose = self.w_lg_nose * self.cg_loc_fus()[1]
            mcg_sys = self.w_sys * self.cg_loc_fus()[2]
            mcg_duct = self.w_duct * self.cg_loc_fus()[3]
            mcg_pax = self.w_pax * self.cg_loc_fus()[4]

            mcg_vect = [mcg_fus, mcg_nose, mcg_sys, mcg_duct, mcg_pax]
            sum_mcg = sum(mcg_vect)

            cg_fuse_group = sum_mcg / sum_w

            components = ["Fuselage", "Nose", "Systems", "Duct", "Pax"]
            file_path = r"C:\Users\tomva\pythonProject\DUUC\data\CG_breakdown"
            os.makedirs(file_path, exist_ok=True)
            input_file_path = os.path.join(file_path, f"DUUC_fuselage_group_cg_breakdown.txt")
            if os.path.exists(input_file_path):
                os.remove(input_file_path)
            with open(input_file_path, "w") as file:
                file.write(print_cg_mass(components, mass_vect, self.cg_loc_fus(), mcg_vect, sum_w, sum_mcg, cg_fuse_group,
                           "fuselage", "DUUC   "))
                file.close()
            return cg_fuse_group, sum_w
        else:
            return None

    def x_cg(self):
        if self.aircraft_type == "conventional":
            x_lemac_guess = 11.5
        elif self.aircraft_type == "DUUC":
            x_lemac_guess = self.x_wing
        else:
            raise ValueError("Unknown aircraft type")
        x_wg_lemac = self.cg_wing_group()[0] - x_lemac_guess
        x_cg_lemac = 0.25 * ref.c_mac_w

        x_lemac = ((self.cg_fuselage_group()[0] - x_cg_lemac) + (self.cg_wing_group()[1] / self.cg_fuselage_group()[1])
                   * (x_wg_lemac - x_cg_lemac))
        #print(f"x_lemac: {x_lemac}")
        x_cg = x_cg_lemac + x_lemac

        x_cg2 = (((self.cg_fuselage_group()[0] * self.cg_fuselage_group()[1]) + (self.cg_wing_group()[0]
                                                                                * self.cg_wing_group()[1])) /
                 (self.cg_wing_group()[1] + self.cg_fuselage_group()[1]))
        return x_cg, x_lemac

    def x_cg_lim(self):
        k = 0.27

        x_for = self.x_cg()[0] - 0.5 * k * self.c_mac_wing
        x_aft = self.x_cg()[0] + 0.5 * k * self.c_mac_wing
        return x_for, x_aft

    """ --------------------------------- CENTER OF GRAVITY Z-DIRECTION ------------------------------------------- """
    def z_loc_components(self):
        z_cg_fus = ref.diameter_fuselage / 2
        z_cg_wing = ref.diameter_fuselage * 0.90
        z_cg_lg = ref.diameter_fuselage * 0.05
        z_cg_pax = z_cg_fus
        z_cg_sys = z_cg_fus
        z_cg_fuel = z_cg_wing
        if self.aircraft_type == "conventional":
            z_cg_ht = ref.diameter_fuselage + ref.b_v
            z_cg_vt = ref.diameter_fuselage + 0.5 * ref.b_v
            z_cg_eng = z_cg_wing
            return [z_cg_fus, z_cg_wing, z_cg_lg, z_cg_pax, z_cg_sys, z_cg_fuel, z_cg_ht, z_cg_vt, z_cg_eng]
        if self.aircraft_type == "DUUC":
            z_cg_pe = self.z_PE
            return [z_cg_fus, z_cg_wing, z_cg_lg, z_cg_pax, z_cg_sys, z_cg_fuel, z_cg_pe]
        else:
            raise ValueError("Wrong aircraft type defined")

    def z_cg(self):
        mz_cg_fus = self.z_loc_components()[0] * self.w_fuse
        mz_cg_wing = self.z_loc_components()[1] * self.w_wing
        mz_cg_lg = self.z_loc_components()[2] * (self.w_lg_nose + self.w_lg_main)
        mz_cg_pax = self.z_loc_components()[3] * self.w_pax
        mz_cg_sys = self.z_loc_components()[4] * self.w_sys
        mz_cg_fuel = self.z_loc_components()[5] * self.w_fuel

        sum_mcg_pre = mz_cg_fus + mz_cg_wing + mz_cg_lg + mz_cg_pax + mz_cg_sys + mz_cg_fuel
        if self.aircraft_type == "conventional":
            mz_cg_ht = self.z_loc_components()[6] * self.w_ht
            mz_cg_vt = self.z_loc_components()[7] * self.w_vt
            mz_cg_eng = self.z_loc_components()[8] * (self.w_eng_wing + self.w_nac_wing) * 2

            m_sum = (self.w_fuse + self.w_wing + self.w_lg_nose + self.w_lg_main + self.w_pax + self.w_sys +
                     self.w_fuel + self.w_ht + self.w_vt + 2 * (self.w_eng_wing + self.w_nac_wing))

            sum_mz_cg = sum_mcg_pre + mz_cg_ht + mz_cg_vt + mz_cg_eng

            z_cg = sum_mz_cg / m_sum
            return z_cg

        if self.aircraft_type == "DUUC":
            mz_cg_pe = self.z_loc_components()[6] * self.w_duct

            m_sum = (self.w_fuse + self.w_wing + self.w_lg_nose + self.w_lg_main + self.w_pax + self.w_sys +
                     self.w_fuel + self.w_duct)

            sum_mz_cg = sum_mcg_pre + mz_cg_pe

            z_cg = sum_mz_cg / m_sum
            return z_cg

        else:
            raise ValueError("Wrong aircraft type defined")


""" Test section"""

if __name__ == "__main__":
    l_fuselage = 21
    x = np.linspace(0, l_fuselage, 100)
    x_wing_group = []
    x_fuselage_group = []
    x_cg = []
    x_cg2 = []
    x_cg3 = []

    x_cg_atr = 11.5586
    x_wg_atr = 11.625
    x_fg_atr = 11.392
    x_wing = 11.2
    x_ac_wing = x_wing + 0.25 * ref.c_mac_w

    for i in range(len(x)):
        cg = CenterOfGravity(w_lg_nose=171,
                             w_lg_main=787,
                             w_eng_wing=500,
                             w_wing=3500,
                             w_fuse=3373,
                             w_nac_wing=250,
                             w_vt=178,
                             w_ht=124,
                             w_duct=3500,
                             w_sys=3113,
                             l_fuselage=l_fuselage,
                             aircraft_type="DUUC",
                             c_mac_wing=ref.c_mac_w,
                             x_wing=x_wing,
                             x_duct=x[i],
                             z_PE=0,
                             w_fuel=config.w_fuel_full_end,
                             w_pax=1000
                             )

        x_wing_group.append(cg.cg_wing_group()[0])
        x_fuselage_group.append(cg.cg_fuselage_group()[0])
        x_cg.append(cg.x_cg()[0])

    for i in range(len(x)):
        cg = CenterOfGravity(w_lg_nose=171,
                             w_lg_main=787,
                             w_eng_wing=500,
                             w_wing=3500,
                             w_fuse=3373,
                             w_nac_wing=250,
                             w_vt=178,
                             w_ht=124,
                             w_duct=3500 - 500,
                             w_sys=3113,
                             l_fuselage=l_fuselage,
                             aircraft_type="DUUC",
                             c_mac_wing=ref.c_mac_w,
                             x_wing=x_wing,
                             x_duct=x[i],
                             z_PE=0,
                             w_fuel=config.w_fuel_full_end,
                             w_pax=1000)

        x_cg2.append(cg.x_cg()[0])

    for i in range(len(x)):
        cg = CenterOfGravity(w_lg_nose=171,
                             w_lg_main=787,
                             w_eng_wing=500,
                             w_wing=3500,
                             w_fuse=3373,
                             w_nac_wing=250,
                             w_vt=178,
                             w_ht=124,
                             w_duct=3500 + 500,
                             w_sys=3113,
                             l_fuselage=l_fuselage,
                             aircraft_type="DUUC",
                             c_mac_wing=ref.c_mac_w,
                             x_wing=x_wing,
                             x_duct=x[i],
                             z_PE=0,
                             w_fuel=config.w_fuel_full_end,
                             w_pax=1000)

        x_cg3.append(cg.x_cg()[0])
    print(x_cg)
    # Create single plot with specified size
    fig, ax0 = plt.subplots(figsize=(12, 4))

    # Load background image
    image = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\DUUC_side_view_tailless.png")

    # Show background image with normalized axes
    ax0.imshow(image, extent=[0, 1, 0, 1], aspect='auto')

    # Plot normalized CG trajectories
    ax0.plot(np.array(x_cg) / l_fuselage, np.array(x) / l_fuselage, label=r"$x_{cg-DUUC}$", color="tab:blue")
    ax0.plot(np.array(x_cg2) / l_fuselage, np.array(x) / l_fuselage, label=r"$x_{cg-DUUC}$: PE - 500 kg",
             color="tab:purple")
    ax0.plot(np.array(x_cg3) / l_fuselage, np.array(x) / l_fuselage, label=r"$x_{cg-DUUC}$: PE + 500 kg",
             color="tab:orange")

    # Span and markers for normalized x_cg range
    ax0.axvline(min(x_cg) / l_fuselage, color="tab:blue", linestyle="dashed")
    ax0.axvline(max(x_cg) / l_fuselage, color="tab:blue", linestyle="dashed")
    ax0.axvspan(min(x_cg) / l_fuselage, max(x_cg) / l_fuselage, color="tab:blue", alpha=0.2)

    # Reference lines
    ax0.axvline(x_cg_atr / l_fuselage, color="tab:green", linestyle="dashed", label=r"$x_{cg-ATR}$")
    ax0.axvline(x_ac_wing / l_fuselage, color="tab:red", label=r"$x_{ac-wing}$")

    # Dotted guideline paths
    ax0.plot([min(x_cg) / l_fuselage, 1], [0, 0], color="black", linestyle="dashed")
    ax0.plot(min(x_cg) / l_fuselage, 0, color="tab:blue", marker="o")
    ax0.plot([np.mean(x_cg) / l_fuselage, 1], [0.5, 0.5], color="black", linestyle="dashed")
    ax0.plot(np.mean(x_cg) / l_fuselage, 0.5, color="tab:blue", marker="o")
    ax0.plot([max(x_cg) / l_fuselage, 1], [1, 1], color="black", linestyle="dashed")
    ax0.plot(max(x_cg) / l_fuselage, 1, color="tab:blue", marker="o")

    # Title and labels
    ax0.set_title("Center of Gravity shift with position of PE")
    ax0.set_xlabel(r"Normalized fuselage location ($x / l_{fuselage}$)")
    ax0.set_ylabel(r"Normalized PE position ($x / l_{fuselage}$)")
    ax0.set_xlim([0, 1])
    ax0.set_ylim([-0.05, 1.05])
    ax0.grid(True)
    ax0.legend()

    plt.tight_layout()
    plt.show()


"""
    l_fuselage = 27
    x = np.linspace(0, l_fuselage, 100)
    x_wing_group = []
    x_fuselage_group = []
    x_cg = []
    x_cg2 = []
    x_cg3 = []
    x_cg4 = []
    x_cg5 = []

    x_cg_atr = 11.5586
    x_wg_atr = 11.625
    x_fg_atr = 11.392
    x_wing = 11.2
    x_ac_wing = x_wing + 0.25 * ref.c_mac_w

    for i in range(len(x)):
        cg = CenterOfGravity(w_lg_nose=171,
                             w_lg_main=787,
                             w_eng_wing=500,
                             w_wing=3500,
                             w_fuse=3373,
                             w_nac_wing=250,
                             w_vt=178,
                             w_ht=124,
                             w_duct=3500,
                             w_sys=3113,
                             l_fuselage=l_fuselage,
                             aircraft_type="DUUC",
                             c_mac_wing=ref.c_mac_w,
                             x_wing=x_wing,
                             x_duct=x[i],
                             z_PE=0,
                             w_fuel=config.w_fuel_full_end,
                             w_pax=1000
                             )

        x_wing_group.append(cg.cg_wing_group()[0])
        x_fuselage_group.append(cg.cg_fuselage_group()[0])
        x_cg.append(cg.x_cg()[0])

    for i in range(len(x)):
        cg = CenterOfGravity(w_lg_nose=171,
                             w_lg_main=787,
                             w_eng_wing=500,
                             w_wing=3500,
                             w_fuse=3373,
                             w_nac_wing=250,
                             w_vt=178,
                             w_ht=124,
                             w_duct=3500 - 500,
                             w_sys=3113,
                             l_fuselage=l_fuselage,
                             aircraft_type="DUUC",
                             c_mac_wing=ref.c_mac_w,
                             x_wing=x_wing,
                             x_duct=x[i],
                             z_PE=0,
                             w_fuel=5000,
                             w_pax=7500)

        x_cg2.append(cg.x_cg()[0])

    for i in range(len(x)):
        cg = CenterOfGravity(w_lg_nose=171,
                             w_lg_main=787,
                             w_eng_wing=500,
                             w_wing=3500,
                             w_fuse=3373,
                             w_nac_wing=250,
                             w_vt=178,
                             w_ht=124,
                             w_duct=3500 + 500,
                             w_sys=3113,
                             l_fuselage=l_fuselage,
                             aircraft_type="DUUC",
                             c_mac_wing=ref.c_mac_w,
                             x_wing=x_wing,
                             x_duct=x[i],
                             z_PE=0,
                             w_fuel=300,
                             w_pax=7500)

        x_cg3.append(cg.x_cg()[0])

    for i in range(len(x)):
        cg = CenterOfGravity(w_lg_nose=171,
                             w_lg_main=787,
                             w_eng_wing=500,
                             w_wing=3500,
                             w_fuse=3373,
                             w_nac_wing=250,
                             w_vt=178,
                             w_ht=124,
                             w_duct=3500 + 500,
                             w_sys=3113,
                             l_fuselage=l_fuselage,
                             aircraft_type="DUUC",
                             c_mac_wing=ref.c_mac_w,
                             x_wing=x_wing,
                             x_duct=x[i],
                             z_PE=0,
                             w_fuel=5000,
                             w_pax=0)

        x_cg4.append(cg.x_cg()[0])

    for i in range(len(x)):
        cg = CenterOfGravity(w_lg_nose=171,
                             w_lg_main=787,
                             w_eng_wing=500,
                             w_wing=3500,
                             w_fuse=3373,
                             w_nac_wing=250,
                             w_vt=178,
                             w_ht=124,
                             w_duct=3500 + 500,
                             w_sys=3113,
                             l_fuselage=l_fuselage,
                             aircraft_type="DUUC",
                             c_mac_wing=ref.c_mac_w,
                             x_wing=x_wing,
                             x_duct=x[i],
                             z_PE=0,
                             w_fuel=300,
                             w_pax=0)

        x_cg5.append(cg.x_cg()[0])

    # Create single plot with specified size
    fig, ax0 = plt.subplots(figsize=(12, 4))

    # Load background imag
"""
#image = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\DUUC_side_view_tailless.png")
"""
    # Show background image with normalized axes
    ax0.imshow(image, extent=[0, 1, 0, 1], aspect='auto')

    # Plot normalized CG trajectories
    ax0.plot(np.array(x_cg) / l_fuselage, np.array(x) / l_fuselage, label=r"OEM")
    ax0.plot(np.array(x_cg2) / l_fuselage, np.array(x) / l_fuselage, label=r"Max payload + max fuel",
             )
    ax0.plot(np.array(x_cg3) / l_fuselage, np.array(x) / l_fuselage, label=r"Max payload + min fuel",
             )
    ax0.plot(np.array(x_cg4) / l_fuselage, np.array(x) / l_fuselage, label=r"No payload + max fuel",
             )
    ax0.plot(np.array(x_cg5) / l_fuselage, np.array(x) / l_fuselage, label=r"No payload + min fuel",
             )

    # Span and markers for normalized x_cg range
    ax0.axvline(min(x_cg) / l_fuselage, color="tab:blue", linestyle="dashed")
    ax0.axvline(max(x_cg) / l_fuselage, color="tab:blue", linestyle="dashed")
    ax0.axvspan(min(x_cg) / l_fuselage, max(x_cg) / l_fuselage, color="tab:blue", alpha=0.2)

    # Reference lines
    #ax0.axvline(x_cg_atr / l_fuselage, color="tab:green", linestyle="dashed", label=r"$x_{cg-ATR}$")
    #ax0.axvline(x_ac_wing / l_fuselage, color="tab:red", label=r"$x_{ac-wing}$")

    # Dotted guideline paths
    #ax0.plot([min(x_cg) / l_fuselage, 1], [0, 0], color="black", linestyle="dashed")
    #ax0.plot(min(x_cg) / l_fuselage, 0, color="tab:blue", marker="o")
    #ax0.plot([np.mean(x_cg) / l_fuselage, 1], [0.5, 0.5], color="black", linestyle="dashed")
    #ax0.plot(np.mean(x_cg) / l_fuselage, 0.5, color="tab:blue", marker="o")
    #ax0.plot([max(x_cg) / l_fuselage, 1], [1, 1], color="black", linestyle="dashed")
    #ax0.plot(max(x_cg) / l_fuselage, 1, color="tab:blue", marker="o")

    # Title and labels
    ax0.set_title("Center of Gravity shift with position of PE")
    ax0.set_xlabel(r"Normalized fuselage location ($x / l_{fuselage}$)")
    ax0.set_ylabel(r"Normalized PE position ($x / l_{fuselage}$)")
    ax0.set_xlim([0, 1])
    ax0.set_ylim([-0.05, 1.05])
    ax0.grid(True)
    ax0.legend()

    plt.tight_layout()
    plt.show()
"""