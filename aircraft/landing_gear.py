import data.atr_reference as ref
import numpy as np


class LandingGear:
    def __init__(self, aircraft_type: str):
        self.aircraft_type = aircraft_type

    """ ----------------------------------- Determine weight Landing Gear ---------------------------------------- """
    def weight(self):
        """ Landing gear mass estimated based on Torenbeek -> output in kg"""
        k = 1.08  # 1.06 for high wing aircraft
        g = 9.81  # gravitational constant [m/s^2]
        w_nose = k * (ref.a_nose + ref.b_nose * ref.MTOW ** 0.75 + ref.d_nose
                      * ref.MTOW ** 1.5)

        w_main = k * (ref.a_main + ref.b_main * ref.MTOW ** 0.75 + ref.d_main
                      * ref.MTOW ** 1.5 + ref.c_main * ref.MTOW)

        return np.round(w_nose, 0), np.round(w_main, 0), np.round(w_nose + w_main, 0)


""" Test section"""
"""
if __name__ == "__main__":
    lg = LandingGear("DUUC")
    print(lg.weight()) """
