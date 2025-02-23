import numpy as np
import flow_conditions


class Wing:
    def __init__(self, b_wing: float, sweep_wing: float, wing_airfoil: str,
                 taper_ratio_wing: float, c_root_wing: float):
        super().__init__()
        self.b_wing = b_wing
        self.sweep_wing = sweep_wing
        self.wing_airfoil = wing_airfoil
        self.tr_wing = taper_ratio_wing
        self.c_root = c_root_wing


