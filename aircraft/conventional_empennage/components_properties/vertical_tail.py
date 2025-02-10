import numpy as np
import constants
import flow_conditions


class VerticalTail:
    def __init__(self, vt_span: float, vt_chord: float, vt_profile: str,
                 vt_taper: float, vt_sweep: float, vt_croot: float, tail_type: str):
        super().__init__()
        self.vt_span = vt_span
        self.vt_chord = vt_chord
        self.vt_profile = vt_profile
        self.tail_type = tail_type
        self.vt_taper = vt_taper
        self.vt_sweep = vt_sweep
        self.vt_croot = vt_croot

    """" Calcultate geometric properties of the horizontal tail"""

    def tip_chord(self):
        c_tip = self.vt_croot * self.vt_taper
        return c_tip

    def area(self):
        """ Area from the vertical stabilizer (based on trapezoid area)"""

        s_vt = 0.5 * self.vt_span * 0.5 * (self.vt_croot + self.vt_chord)
        print(f'Area of the vertical stabilizer = {2* s_vt} [m^2]')
        return s_vt

    def x_vt_tail(self):
        """ Taking the leading edge of the chord profile as a reference"""
        sweep_rad = np.radians(self.vt_sweep)
        x_tip_tail = 0.5 * self.vt_span * np.sin(sweep_rad)
        return x_tip_tail

    def z_vt_tail(self):
        """ depending on the tail type determine the z-location of the tail"""
        if self.tail_type == 'T':
            z_tail = self.vt_span
            return z_tail

        else:
            z_tail = 0
            return z_tail
