import constants


def reynolds(rho, u, char_len):
    re = (rho * u * char_len) / constants.mu_air
    return re
