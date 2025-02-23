import numpy as np
import os

""" Input section"""
prop_airfoil = "ARAD8"


def ct_function(a, glauert=False):
    ct = np.zeros(np.shape(a))
    ct = 4 * a * (1 - a)
    if glauert:
        ct1 = 1.816
        a1 = 1 - np.sqrt(ct1) / 2
        ct[a > a1] = ct1 - 4 * (np.sqrt(ct1) - 1) * (1 - a[a > a1])
    return ct


def a_induction(CT):
    """ This function calculates the induction factor 'a' as a function of thrust
    coefficient CT  including Glauert correction """
    a = np.zeros(np.shape(CT))
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    a[CT >= CT2] = 1 + (CT[CT >= CT2] - CT1) / (4 * (np.sqrt(CT1) - 1))
    a[CT < CT2] = 0.5 - 0.5 * np.sqrt(1 - CT[CT < CT2])
    return a

# plot CT as a function of induction "a", with and without Glauert correction
# define a as a range
# a = np.arange(-.5,1,.01)
# CTmom = ct_function(a) # CT without correction
# CTglauert = ct_function(a, True) # CT with Glauert's correction
# a2 = a_induction(CTglauert)

# fig1 = plt.figure(figsize=(12, 6))
# plt.plot(a, CTmom, 'k-', label='$C_T$')
# plt.plot(a, CTglauert, 'b--', label='$C_T$ Glauert')
# plt.plot(a, CTglauert*(1-a), 'g--', label='$C_P$ Glauert')
# plt.xlabel('a')
# plt.ylabel(r'$C_T$ and $C_P$')
# plt.grid()
# plt.legend()


def prandtl_correction(r_R, rootradius_r, tipradius_r, tsr, n_blades, axial_induction):
    """ This function calculates the combined tip and root Prandtl correction at a given radial
    position 'r_R' (non-dimensioned by rotor radius), given a root and tip radius
    (also non-dimensioned), a tip speed ratio TSR, the number lf blades NBlades and the axial
    induction factor """

    temp1 = -n_blades/2*(tipradius_r-r_R)/r_R*np.sqrt(1 + ((tsr*r_R)**2)/((1-axial_induction)**2))
    F_tip = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    F_tip[np.isnan(F_tip)] = 0
    temp1 = n_blades/2*(rootradius_r-r_R)/r_R*np.sqrt(1 + ((tsr*r_R)**2)/((1-axial_induction)**2))
    F_root = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    F_root[np.isnan(F_root)] = 0
    return F_root*F_tip, F_tip, F_root


# plot Prandtl tip, root and combined correction for a number of blades and induction 'a', over the non-dimensioned radius
# r_R = np.arange(0.1, 1, .01)
# a = np.zeros(np.shape(r_R))+0.3
# Prandtl, Prandtltip, Prandtlroot = prandtl_correction(r_R, 0.1, 1, 7, 3, a)

# fig1 = plt.figure(figsize=(12, 6))
# plt.plot(r_R, Prandtl, 'r-', label='Prandtl')
# plt.plot(r_R, Prandtltip, 'g.', label='Prandtl tip')
# plt.plot(r_R, Prandtlroot, 'b.', label='Prandtl root')
# plt.xlabel('r/R')
# plt.legend()


def data_polar(prop_airfoil):
    filepath = os.path.join(r"C:\Users\tomva\pythonProject\DUUC\data\Polars",
                            f"{prop_airfoil}_polar.txt")
    # Normalize the filepath to ensure it resolves correctly
    filepath = os.path.abspath(filepath)

    data = np.loadtxt(filepath, skiprows=2)
    polar_alpha = data[:, 0]
    polar_cl = data[:, 1]
    polar_cd = data[:, 2]
    return polar_alpha, polar_cl, polar_cd


# fig, axs = plt.subplots(1, 2, figsize=(12, 6))
# axs[0].plot(polar[0], polar[1])
# axs[0].set_xlim([-30, 30])
# axs[0].set_xlabel(r'$\alpha$')
# axs[0].set_ylabel(r'$C_l$')
# axs[0].grid()
# axs[1].plot(polar[2], polar[1])
# axs[1].set_xlim([0, .1])
# axs[1].set_xlabel(r'$C_d$')
# axs[1].grid()


# define function to determine load in the blade element
def load_blade_element(v_norm, v_tan, r_R, chord, twist, polar_alpha, polar_cl, polar_cd):
    """
    calculates the load in the blade element
    """
    v_mag2 = np.clip(v_norm ** 2 + v_tan ** 2, 0, 1e6)
    inflow_angle = np.arctan2(v_norm, v_tan)
    print(f"inflow angle: {inflow_angle}")
    alpha = twist - np.rad2deg(inflow_angle)
    print(f"alpha: {alpha}")
    cl = np.interp(alpha, polar_alpha, polar_cl)
    cd = np.interp(alpha, polar_alpha, polar_cd)
    lift = 0.5 * v_mag2 * cl * chord
    drag = 0.5 * v_mag2 * cd * chord
    f_norm = lift * np.cos(inflow_angle) + drag * np.sin(inflow_angle)
    f_tan = lift * np.sin(inflow_angle) - drag * np.cos(inflow_angle)
    gamma = 0.5 * np.sqrt(v_mag2) * cl * chord
    # print(f_norm, f_tan, gamma)
    return f_norm, f_tan, gamma


def solve_streamtube(u_inf, r1_r, r2_r, rootradius_r, tipradius_r, omega, radius, n_blades, chord,
                     twist, polar_alpha, polar_cl, polar_cd):
    """
    solve balance of momentum between blade element load and loading in the stream tube
    input variables:
    U_inf - wind speed at infinity
    r1_R,r2_R - edges of blade element, in fraction of Radius ;
    root radius_R, tip radius_R - location of blade root and tip, in fraction of Radius ;
    Radius is the rotor radius
    Omega -rotational velocity
    NBlades - number of blades in rotor
    """
    area = np.pi * ((r2_r * radius) ** 2 - (r1_r * radius) ** 2)  # area stream tube
    r_R = (r1_r + r2_r) / 2  # centroide

    # initialize variables
    a = 0.0  # axial induction
    aline = 0.0  # tangential induction factor

    n_iterations = 100
    error_iterations = 0.00001  # error limit for iteration process, in absolute value of induction

    for i in range(n_iterations):
        # ///////////////////////////////////////////////////////////////////////
        # // this is the block "Calculate velocity and loads at blade element"
        # ///////////////////////////////////////////////////////////////////////
        u_rotor = u_inf * (1 - a)  # axial velocity at rotor
        u_tan = (1 + aline) * omega * r_R * radius  # tangential velocity at rotor
        # calculate loads in blade segment in 2D (N/m)
        fnorm, ftan, gamma = load_blade_element(u_rotor, u_tan, r_R, chord, twist, polar_alpha,
                                                polar_cl, polar_cd)
        load3Daxial = fnorm * radius * (r2_r - r1_r) * n_blades  # 3D force in axial direction
        # load3D tan =loads[1]*Radius*(r2_R-r1_R)*NBlades # 3D force in azimuthal/
        # tangential direction (not used here)

        # ///////////////////////////////////////////////////////////////////////
        # //the block "Calculate velocity and loads at blade element" is done
        # ///////////////////////////////////////////////////////////////////////

        # ///////////////////////////////////////////////////////////////////////
        # // this is the block "Calculate new estimate of axial and azimuthal induction"
        # ///////////////////////////////////////////////////////////////////////
        # // calculate thrust coefficient at the stream tube
        CT = load3Daxial / (0.5 * area * u_inf ** 2)

        # calculate new axial induction, accounting for Glauert's correction
        anew = a_induction(CT)

        # correct new axial induction with Prandtl's correction
        Prandtl, Prandtltip, Prandtlroot = prandtl_correction(r_R, rootradius_r, tipradius_r,
                                                              omega * radius / u_inf, n_blades,
                                                              anew)
        if (Prandtl < 0.0001):
            Prandtl = 0.0001  # avoid divide by zero
        anew = anew / Prandtl  # correct estimate of axial induction
        a = 0.75 * a + 0.25 * anew  # for improving convergence, weigh current and
        # previous iteration of axial induction

        # calculate axi mu induction
        aline = ftan * n_blades / (2 * np.pi * u_inf * (1 - a) * omega * 2 * (r_R * radius) ** 2)
        aline = aline / Prandtl  # correct estimate of azimuthal induction with Prandtl's correction
        # ///////////////////////////////////////////////////////////////////////////
        # // end of the block "Calculate new estimate of axial and azimuthal induction"
        # ///////////////////////////////////////////////////////////////////////

        # // test convergence of solution, by checking convergence of axial induction
        if (np.abs(a - anew) < error_iterations):
            # print("iterations")
            # print(i)
            break

    return [a, aline, r_R, fnorm, ftan, gamma]


def bem_run(pitch, speed, radius, n_blades, prop_airfoil, rpm):
    # define the blade geometry
    delta_r_R = .01
    r_R = np.arange(0.2, 1 + delta_r_R / 2, delta_r_R)

    # blade shape
    chord_distribution = 3 * (1 - r_R) + 1  # meters
    twist_distribution = -14 * (1 - r_R) + pitch  # degrees

    # define flow conditions
    omega = rpm * 2 * np.pi / 60
    tsr = omega * radius / speed

    tip_r = 1
    root_r = 0.2

    polar = data_polar(prop_airfoil)

    # solve BEM model
    results = np.zeros([len(r_R) - 1, 6])

    for i in range(len(r_R) - 1):
        chord = np.interp((r_R[i] + r_R[i + 1]) / 2, r_R, chord_distribution)
        twist = np.interp((r_R[i] + r_R[i + 1]) / 2, r_R, twist_distribution)

        results[i, :] = solve_streamtube(speed, r_R[i], r_R[i + 1], root_r, tip_r, omega,
                                         radius, n_blades, chord, twist, polar[0],
                                         polar[1], polar[2])

    areas = (r_R[1:]**2-r_R[:-1]**2)*np.pi*radius**2
    dr = (r_R[1:]-r_R[:-1])*radius
    CT = np.sum(dr*results[:, 3]*n_blades/(0.5*speed**2*np.pi*radius**2))
    CP = np.sum(dr*results[:, 4]*results[:, 2]*n_blades*radius * omega/(0.5*speed**3*np.pi*radius))

    # Calculate total thrust and torque
    total_thrust = np.sum(dr * results[:, 3] * n_blades)
    total_torque = np.sum(dr * results[:, 4] * results[:, 2] * n_blades * radius)

    print("\n----- BEM calculator -----")
    print(f"CT = {np.round(CT, 2)} [-]")
    print(f"CP = {np.round(CP, 2)} [-]")
    print(f"Total thrust = {np.round(total_thrust)} [N]")
    print(f"Total torque = {np.round(total_torque)} [Nm]")
    return CT, CP, total_thrust, total_torque


a = bem_run(46, 148, 1.8, 3, "ARAD8", 5000)
