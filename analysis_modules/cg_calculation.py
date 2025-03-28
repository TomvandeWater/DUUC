from analysis_modules.plotting_functions import print_cg_mass
import data.atr_reference as ref


class CenterOfGravity:
    """ determine center of gravity based on method of Scholz 2008, note: input of weight in [kg] and distance in [m]"""
    def __init__(self, w_lg_nose: float, w_lg_main: float, w_eng_wing: float, w_wing: float, w_fuse: float,
                 w_nac_wing: float, w_vt: float, w_ht: float, w_duct: float, w_sys: float, l_fuselage: float,
                 aircraft_type: str, c_mac_wing: float):
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
            cg_wing = 12.24
            cg_main = 12.429
            return [cg_wing, cg_main]
        else:
            return None

    def cg_loc_fus(self):
        if self.aircraft_type == "conventional":
            cg_fus = 0.39 * self.l_fuselage
            cg_vt = 23.99
            cg_ht = 23.99
            cg_nose = 1.75
            cg_sys = 11.5
            return [cg_fus, cg_vt, cg_ht, cg_nose, cg_sys]
        if self.aircraft_type == "DUUC":
            cg_fus = 0.5 * self.l_fuselage
            cg_nose = 1.75
            cg_sys = 11.5
            cg_duct = 23.99  # will be variable later
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
            # print_cg_mass(components, mass_vect, self.cg_loc_wing(), mcg_vect, sum_w, sum_mcg, cg_wing_group,
                          # "wing    ", self.aircraft_type)

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
            # print_cg_mass(components, mass_vect, self.cg_loc_wing(), mcg_vect, sum_w, sum_mcg, cg_wing_group,
                          # "wing     ", self.aircraft_type)

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
            # print_cg_mass(components, mass_vect, self.cg_loc_fus(), mcg_vect, sum_w, sum_mcg, cg_fuse_group,
                          # "fuselage", self.aircraft_type)

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
            # print_cg_mass(components, mass_vect, self.cg_loc_fus(), mcg_vect, sum_w, sum_mcg, cg_fuse_group,
                          # "fuselage", self.aircraft_type)

            return cg_fuse_group, sum_w
        else:
            return None

    def x_cg(self):
        x_lemac_guess = 11.5
        x_wg_lemac = self.cg_wing_group()[0] - x_lemac_guess
        x_cg_lemac = 0.25 * ref.c_mac_w

        x_lemac = (self.cg_fuselage_group()[0] - x_cg_lemac + (self.cg_wing_group()[1] / self.cg_fuselage_group()[1])
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
    cg = CenterOfGravity(w_lg_nose=171,
                         w_lg_main=700,
                         w_eng_wing=250,
                         w_wing=3500,
                         w_fuse=3000,
                         w_nac_wing=200,
                         w_vt=1002,
                         w_ht=1001,
                         w_duct=2000,
                         w_sys=1000,
                         l_fuselage=ref.l_cab + ref.l_tail + ref.l_cockpit,
                         aircraft_type="DUUC")
    print(cg.cg_wing_group())
    print(cg.cg_fuselage_group())
