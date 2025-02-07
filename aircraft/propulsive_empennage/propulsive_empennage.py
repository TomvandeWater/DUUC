from aircraft.propulsive_empennage.component_properties.pylon import Pylon


class PropulsiveEmpennage:
    def __init__(self, pylon_length: float, pylon_chord: float,
                 pylon_profile: str, cant_angle: float, alpha: float):

        self.pylon = Pylon(pylon_length, pylon_chord, pylon_profile,
                           cant_angle, alpha)

    def thrust_pe(self):
        thrust = 150
        print(f"\nThrust of the pylon is", thrust)
        return thrust
