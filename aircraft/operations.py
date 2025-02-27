import config
import data.atr_reference as ref


class Operations:
    def __init__(self, pax: float, MTOM: float):
        self.pax = pax
        self.MTOW = MTOM

    """ ------------------------------------- Determine weights ------------------------------------------------- """
    def weight_pax(self):
        g = 9.81  # gravitational constant

        w_pass = self.pax * config.w_pax
        w_add = self.pax * 15  # additional 15 kg per passenger for seats etc.

        w_pax = (w_pass + w_add) * g
        return w_pax

    def weight_sys(self):
        """ considers additional weight of systems -> output in newtons"""
        k_equip = 0.11  # from text nita
        g = 9.81  # gravitational constant

        w_sys = k_equip * self.MTOW * g  # weight of the control actuation ignored for now.
        return w_sys

    @staticmethod
    def weight_fuel():
        g = 9.81  # gravitational constant

        w_fuel = 1989 * g  # reference fuel mass
        return w_fuel


""" Test section """

if __name__ == "__main__":
    operations = Operations(config.n_pax, ref.MTOW)
    print(f"weight pax: {operations.weight_pax()}")
    print(f"weight sys: {operations.weight_sys()}")
    print(f"weight fuel: {operations.weight_fuel()}")
