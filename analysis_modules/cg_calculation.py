import numpy as np
import matplotlib.pyplot as plt
from analysis_modules.plotting_functions import print_cg_mass
import data.atr_reference as ref
import os
import matplotlib.gridspec as gridspec


class CenterOfGravity:
    """ determine center of gravity based on method of Scholz 2008, note: input of weight in [kg] and distance in [m]"""
    def __init__(self, w_lg_nose: float, w_lg_main: float, w_eng_wing: float, w_wing: float, w_fuse: float,
                 w_nac_wing: float, w_vt: float, w_ht: float, w_duct: float, w_sys: float, l_fuselage: float,
                 aircraft_type: str, c_mac_wing: float, x_duct: float, x_wing: float):
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

    """ ---------------------------- determine center of gravity of each component -------------------------------- """
    def cg_loc_wing(self):
        if self.aircraft_type == "conventional":
            """ Data extracted from reference values ATR72-600"""
            cg_wing = 12.24
            cg_main = 12.429
            cg_eng_wing = 10.12
            cg_nac_wing = 10.12
            return [cg_wing, cg_main, cg_eng_wing, cg_nac_wing]
        if self.aircraft_type == "DUUC":
            cg_wing = self.x_wing + self.c_mac_wing / 2
            cg_main = cg_wing + 0.25
            return [cg_wing, cg_main]
        else:
            return None

    def cg_loc_fus(self):
        if self.aircraft_type == "conventional":
            cg_fus = 0.45 * self.l_fuselage
            cg_vt = 23.99
            cg_ht = 23.99
            cg_nose = 1.75
            cg_sys = 11.5
            return [cg_fus, cg_vt, cg_ht, cg_nose, cg_sys]
        if self.aircraft_type == "DUUC":
            cg_fus = 0.45 * self.l_fuselage
            cg_nose = 1.75
            cg_sys = 11.5
            cg_duct = self.x_duct  # will be variable later
            return [cg_fus, cg_nose, cg_sys, cg_duct]

    def cg_wing_group(self):
        if self.aircraft_type == 'conventional':
            mass_vect = [self.w_wing, self.w_lg_main, self.w_eng_wing * 2, self.w_nac_wing * 2]

            sum_w = sum(mass_vect)

            mcg_wing = self.w_wing * self.cg_loc_wing()[0]
            mcg_main = self.w_lg_main * self.cg_loc_wing()[1]
            mcg_eng_wing = self.w_eng_wing * 2 * self.cg_loc_wing()[2]
            mcg_nac_wing = self.w_nac_wing * 2 * self.cg_loc_wing()[3]

            mcg_vect = [mcg_wing, mcg_main, mcg_eng_wing, mcg_nac_wing]

            sum_mcg = sum(mcg_vect)

            cg_wing_group = sum_mcg / sum_w

            components = ["Wing", "Main", "Engine", "Nacelle"]
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
            mass_vect = [self.w_wing, self.w_lg_main]

            sum_w = sum(mass_vect)

            mcg_wing = self.w_wing * self.cg_loc_wing()[0]
            mcg_main = self.w_lg_main * self.cg_loc_wing()[1]

            mcg_vect = [mcg_wing, mcg_main]

            sum_mcg = sum(mcg_vect)

            cg_wing_group = sum_mcg / sum_w

            components = ["Wing", "Main"]
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
            mass_vect = [self.w_fuse, self.w_vt, self.w_ht, self.w_lg_nose, self.w_sys]

            sum_w = sum(mass_vect)

            mcg_fus = self.w_fuse * self.cg_loc_fus()[0]
            mcg_vt = self.w_vt * self.cg_loc_fus()[1]
            mcg_ht = self.w_ht * self.cg_loc_fus()[2]
            mcg_nose = self.w_lg_nose * self.cg_loc_fus()[3]
            mcg_sys = self.w_sys * self.cg_loc_fus()[4]

            mcg_vect = [mcg_fus, mcg_vt, mcg_ht, mcg_nose, mcg_sys]

            sum_mcg = sum(mcg_vect)

            cg_fuse_group = sum_mcg / sum_w

            components = ["Fuselage", "Vertical Tail", "Horizontal Tail", "Nose", "Systems"]

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
            mass_vect = [self.w_fuse, self.w_lg_nose, self.w_sys, self.w_duct]
            sum_w = sum(mass_vect)

            mcg_fus = self.w_fuse * self.cg_loc_fus()[0]
            mcg_nose = self.w_lg_nose * self.cg_loc_fus()[1]
            mcg_sys = self.w_sys * self.cg_loc_fus()[2]
            mcg_duct = self.w_duct * self.cg_loc_fus()[3]

            mcg_vect = [mcg_fus, mcg_nose, mcg_sys, mcg_duct]
            sum_mcg = sum(mcg_vect)

            cg_fuse_group = sum_mcg / sum_w

            components = ["Fuselage", "Nose", "Systems", "Duct"]
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
        x_lemac_guess = 11.5
        x_wg_lemac = self.cg_wing_group()[0] - x_lemac_guess
        x_cg_lemac = 0.25 * ref.c_mac_w

        x_lemac = ((self.cg_fuselage_group()[0] - x_cg_lemac) + (self.cg_wing_group()[1] / self.cg_fuselage_group()[1])
                   * (x_wg_lemac - x_cg_lemac))

        x_cg = x_cg_lemac + x_lemac
        return x_cg, x_lemac

    def x_cg_lim(self):
        k = 0.27

        x_for = self.x_cg()[0] - 0.5 * k * self.c_mac_wing
        x_aft = self.x_cg()[0] + 0.5 * k * self.c_mac_wing
        return x_for, x_aft


""" Test section"""

if __name__ == "__main__":
    l_fuselage = ref.l_cab + ref.l_tail + ref.l_cockpit
    x = np.linspace(0, l_fuselage, 100)
    x_wing_group = []
    x_fuselage_group = []
    x_cg = []

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
                             w_duct=1750,
                             w_sys=3113,
                             l_fuselage=l_fuselage,
                             aircraft_type="DUUC",
                             c_mac_wing=ref.c_mac_w,
                             x_wing=x_wing,
                             x_duct=x[i])

        x_wing_group.append(cg.cg_wing_group()[0])
        x_fuselage_group.append(cg.cg_fuselage_group()[0])
        x_cg.append(cg.x_cg()[0])

    fig = plt.figure(figsize=(18, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])  # Left is 3x wider than right

    # Left subplot
    ax0 = plt.subplot(gs[0])
    image = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\DUUC_side_view_tailless.png")

    # Plot background image (aircraft side view)
    ax0.imshow(image, extent=[0, l_fuselage, 0, l_fuselage], aspect='auto')
    #ax0.plot(x_wing_group, x, label=r"$x_{wing group}$", color="tab:blue")
    #ax0.plot(x_fuselage_group, x, label=r"$x_{fuselage group}$", color="tab:orange")
    ax0.plot(x_cg, x, label=r"$x_{cg-DUUC}$", color="tab:blue")
    ax0.axvline(min(x_cg), color="tab:blue", linestyle="dashed")
    ax0.axvline(max(x_cg), color="tab:blue", linestyle="dashed")
    ax0.axvspan(min(x_cg), max(x_cg), color="tab:blue", alpha=0.2)

    ax0.axvline(x_cg_atr, color="tab:green", linestyle="dashed", label="$x_{cg-atr}$")
    ax0.axvline(x_ac_wing, color="tab:red", label=r"$x_{ac-wing}$")
    #ax0.axvline(x_wg_atr, color="tab:blue", linestyle="dashed", alpha=0.5, label="$x_{fg-atr}$")
    #ax0.axvline(x_fg_atr, color="tab:orange", linestyle="dashed", alpha=0.5, label="$x_{wg-atr}$")

    ax0.plot([min(x_cg), l_fuselage], [0, 0], color="black", linestyle="dashed")
    ax0.plot(min(x_cg), 0, color="tab:blue", marker="o")
    ax0.plot([np.mean(x_cg), l_fuselage], [0.5 * l_fuselage, 0.5 * l_fuselage], color="black", linestyle="dashed")
    ax0.plot(np.mean(x_cg), 0.5 * l_fuselage, color="tab:blue", marker="o")
    ax0.plot([max(x_cg), l_fuselage], [l_fuselage, l_fuselage], color="black", linestyle="dashed")
    ax0.plot(max(x_cg), l_fuselage, color="tab:blue", marker="o")

    ax0.set_title("Center of Gravity Shift with Position of PE")
    ax0.set_xlabel(r"Fuselage Location [m]")
    ax0.set_ylabel(f"x-location of the PE on the fuselage [m]")
    ax0.grid(True)
    ax0.set_xlim([0, l_fuselage])
    ax0.set_ylim([-1, l_fuselage*1.05])
    ax0.legend()

    image2 = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\DUUC_diff_loc.png")
    # Right subplot - put whatever figure or image you want here
    ax1 = plt.subplot(gs[1])
    # Example: draw a placeholder box or text
    ax1.imshow(image2, aspect='auto')
    ax1.axis('off')  # Hide axes for a cleaner look

    plt.tight_layout()
    plt.show()



