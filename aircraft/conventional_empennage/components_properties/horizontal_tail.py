import numpy as np
import constants
import flow_conditions


class HorizontalTail:
    def __init__(self, ht_span: float, ht_chord: float, ht_profile: str,
                 ht_taper: float, ht_sweep: float, ht_croot: float):
        super().__init__()
        self.ht_span = ht_span
        self.ht_chord = ht_chord
        self.ht_profile = ht_profile
        self.ht_taper = ht_taper
        self.ht_sweep = ht_sweep
        self.ht_croot = ht_croot

    """" Calcultate geometric properties of the horizontal tail"""

    def tip_chord(self):
        c_tip = self.ht_croot * self.ht_taper
        return c_tip

    def area(self):
        """ Area from one horizontal stabilizer (based on trapezoid area)"""

        s_ht = 0.5 * self.ht_span * 0.5 * (self.ht_croot + self.ht_chord)
        print(f'Area of the horizontail stabilizer = {2* s_ht} [m^2]')
        return s_ht

    def x_ht_tail(self):
        """ Taking the leading edge of the chord profile as a reference"""
        sweep_rad = np.radians(self.ht_sweep)
        x_tip_tail = 0.5 * self.ht_span * np.sin(sweep_rad)

        return x_tip_tail
