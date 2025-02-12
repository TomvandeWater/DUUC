from aircraft.conventional_empennage.components_properties.horizontal_tail import HorizontalTail
from aircraft.conventional_empennage.components_properties.vertical_tail import VerticalTail


class ConventionalEmpennage:
    def __init__(self, ht_span: float, ht_chord: float, ht_profile: str, ht_taper: float,
                 ht_sweep: float, ht_croot: float, vt_span: float, vt_chord: float, vt_profile: str,
                 vt_taper: float, vt_sweep: float, vt_croot: float, tail_type: str):

        self.horizontal_tail = HorizontalTail(ht_span, ht_chord, ht_profile, ht_taper,
                                              ht_sweep, ht_croot)

        self.vertical_tail = VerticalTail(vt_span, vt_chord, vt_profile, vt_taper, vt_sweep,
                                          vt_croot, tail_type)

    def fx(self):
        fx_emp = 0
        return fx_emp

    def fy(self):
        fy_emp = 1
        return fy_emp

    def fz(self):
        fz_emp = 0
        return fz_emp