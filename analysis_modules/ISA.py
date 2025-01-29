import constants


def air_density_isa(altitude):
    """
    Calculate the air density based on the International Standard Atmosphere (ISA) model.
    :constants: are imported from constants file
    :param altitude: Altitude in meters
    :return: Air density in kg/m³
    """

    if altitude < 11000:  # Troposphere (up to 11 km)
        T = constants.T0 - constants.L * altitude
        # Temperature at altitude (K)
        P = constants.P0 * (T / constants.T0) ** (constants.g /
                                                  (constants.R * constants.L))
        # Pressure at altitude (Pa)
    else:  # Stratosphere simplification (above 11 km, up to 20 km)
        T = 216.65  # Constant temperature in the lower stratosphere
        P = constants.P0 * 0.22336 * (2.718 ** (-constants.g
                                                * (altitude - 11000) /
                                                (constants.R * T)))

    rho = P / (constants.R * T)  # Air density (kg/m³)
    return rho
