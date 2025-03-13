import numpy as np
import os

data_folder = r"C:\Users\tomva\pythonProject\DUUC\data"
polar_folder = r"C:\Users\tomva\pythonProject\DUUC\data\Polars"
output_folder = r"C:\Users\tomva\pythonProject\DUUC\data"


# Define a function to clear the screen
def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')


# Subroutine for external commands
def system(command):
    os.system(command)


# Main loop
def propeller_analysis(prop_name: str, choice_ducted: str, choice_induced: str):
    # Clear screen
    #clear_screen()

    # Constants
    max_psi_i = 30
    max_rr_i = 25

    # Default values
    new_case = 0
    zero = 1.0e-20
    pi = 4.0 * np.arctan(1.0)
    tempkzero = 288.16  # at zero S.A
    mhuzero = 1.7894e-5  # idem
    gnu_output = 1  # default value for gnu output
    dimxx = 1
    nix = 1
    niy = 1
    use_induced = False
    ducted = False
    scale_cd = 1.0

    # Directory for gnuplot macros
    if new_case == 0:
        gnuplotmacros = '~/surfdrive/5A-UTILITIES-1/fortran/prop1b/gnuplotmacros/'
        # print(gnuplotmacros)
    print('---------- BEM --------')
    print(' *---------------------------------------------*')
    print(' |                 P R O P 1 B                 |')
    print(' |                                             |')
    print(' |                 version 2.2                 |')
    print(' |                                             |')
    print(' |            Author: L.L.M. Veldhuis          |')
    print(' |                                             |')
    print(' |        Delft University of Technology       |')
    print(' |       Faculty of Aerospace Engineering      |')
    print(' |                                             |')
    print(' *---------------------------------------------*')

    # Input file with blade data
    while True:
        blade_data = os.path.join(data_folder, f"{prop_name}_data.txt")
        try:
            with open(blade_data, 'r') as f:
                propid = f.readline().strip()
                f.readline()  # empty line
                diam = float(f.readline().strip())
                nrblades = int(f.readline().strip())
                tempc = float(f.readline().strip())
                rho = float(f.readline().strip())
                alphap = float(f.readline().strip())
                tas = float(f.readline().strip())
                advratio = float(f.readline().strip())
                beta75 = float(f.readline().strip())
                nrb_sections = int(f.readline().strip())

                rR = np.zeros(nrb_sections)
                cr = np.zeros(nrb_sections)
                beta0 = np.zeros(nrb_sections)
                for i in range(nrb_sections):
                    line = f.readline().strip().split()
                    rR[i] = float(line[0])
                    cr[i] = float(line[1])
                    beta0[i] = float(line[2])

                # Tip r over R should be slightly smaller than 1.0
                if abs(rR[nrb_sections - 1] - 1.0) < zero:
                    rR[nrb_sections - 1] = 0.99

                samediff = f.readline().strip()
                if samediff == 'same':
                    airfoilfile_name = f.readline().strip()
                    airfoilfile = os.path.join(polar_folder, f"{prop_name}_polar.txt")

                nr_alfa = np.zeros(nrb_sections, dtype=int)
                ma1_or_re1 = np.zeros(nrb_sections)
                ma2_or_re2 = np.zeros(nrb_sections)
                alfa = np.zeros((nrb_sections, 28))
                cl1 = np.zeros((nrb_sections, 28))
                cd1 = np.zeros((nrb_sections, 28))
                cl2 = np.zeros((nrb_sections, 28))
                cd2 = np.zeros((nrb_sections, 28))
                re1 = np.zeros(nrb_sections)
                re2 = np.zeros(nrb_sections)
                mach1 = np.zeros(nrb_sections)
                mach2 = np.zeros(nrb_sections)

                for i in range(nrb_sections):
                    if samediff == 'same':
                        airfoil_file = airfoilfile
                    else:
                        airfoil_file = f.readline().strip()

                    try:
                        with open(airfoil_file, 'r') as af:
                            nr_alfa[i] = int(af.readline().strip())
                            ma1_or_re1[i] = float(af.readline().strip())
                            for j in range(nr_alfa[i]):
                                line = af.readline().strip().split()
                                alfa[i, j] = float(line[0])
                                cl1[i, j] = float(line[1])
                                cd1[i, j] = float(line[2])
                            ma2_or_re2[i] = float(af.readline().strip())
                            for j in range(nr_alfa[i]):
                                line = af.readline().strip().split()
                                alfa[i, j] = float(line[0])
                                cl2[i, j] = float(line[1])
                                cd2[i, j] = float(line[2])
                    except FileNotFoundError:
                        print(f"Error: Airfoil file {airfoil_file} not found.")
                        proceed = False
                        break  # Exit the section loop
                    except Exception as e:
                        print(f"Error reading airfoil file {airfoil_file}: {e}")
                        proceed = False
                        break

                    re1[i] = -1.
                    re2[i] = -1.
                    mach1[i] = -1.
                    mach2[i] = -1.
                    if (ma1_or_re1[i] > 1 and ma2_or_re2[i] > 1):
                        re1[i] = ma1_or_re1[i]
                        re2[i] = ma2_or_re2[i]
                    elif (ma1_or_re1[i] < 1 and ma2_or_re2[i] < 1):
                        mach1[i] = ma1_or_re1[i]
                        mach2[i] = ma2_or_re2[i]
                    else:
                        print(' WRONG DATA IN DATAFILE')
                        print(' Mach_or_Re -1 and -2 not both Mach or Re.')
                        proceed = False
                        break  # Exit inner loop, go back to file input
                else:
                    break  # Exit the loop if file reading was successful
                continue
            break  # Exit the loop if file reading was successful

        except FileNotFoundError:
            print("Error: File not found. Please enter a valid filename.")
        except Exception as e:
            print(f"An error occurred: {e}")

    # Since J, n and tas are coupled, first calculate rpm based on input
    rpm = (tas / (advratio * diam)) * 60  # RPM is calculated

    # Ducted propeller?
    if new_case == 0:
        pass
    # answer = input(' Ducted propeller ? (y/n) :')
    answer = 'y'
    if answer == 'y':
        ducted = True

    # Use wing induced velocity data?
    use_induced = False
    # answer = input(' Use (wing) induced velocity data ? (y/n) :')
    answer = choice_induced
    if answer == 'y':
        use_induced = True
        file_induced = input(' Enter file containing induced velocities :')

        # Check whether the file contains the complete velocity vector or just
        # the axial induced velocity
        print(' Select file contents:')
        print(' (1) psi,r/R,Vx,Vy,Vz')
        print(' (2) psi,r/R,Vx (Vy=Vz=0)')
        ifilecontents = int(input(' Your choice:'))

        try:
            with open(file_induced, 'r') as f:
                nr_psi_i, nr_rr_i = map(int, f.readline().strip().split())

                if nr_psi_i > max_psi_i or nr_rr_i > max_rr_i:
                    print(' Too many induced data points')
                    exit()

                psi_i = np.zeros(nr_psi_i)  # Corrected size
                rr_i = np.zeros(nr_rr_i)  # Corrected size
                vx_i = np.zeros((nr_psi_i, nr_rr_i))
                vy_i = np.zeros((nr_psi_i, nr_rr_i))
                vz_i = np.zeros((nr_psi_i, nr_rr_i))

                if ifilecontents == 1:
                    for i in range(nr_psi_i):
                        line = f.readline().strip().split()
                        psi_i[i], rr_i[i], vx_i[i, i], vy_i[i, i], vz_i[i, i] = map(float, line)  # Corrected indexing
                else:
                    for i in range(nr_psi_i):
                        line = f.readline().strip().split()
                        psi_i[i], rr_i[i], vx_i[i, i] = map(float, line)  # Corrected indexing
                        vy_i[i, i] = 0.0
                        vz_i[i, i] = 0.0
        except FileNotFoundError:
            print("Error: Induced velocity file not found.")
        except Exception as e:
            print(f"An error occurred while reading induced velocity data: {e}")

    # Interactive input data
    while True:
        #clear_screen()

        print('  Input data:')
        print(f'  1) Diameter of prop (m)        = {diam:10.3f}')
        print(f'  2) Beta at 75% (degr.)         = {beta75:9.2f}')
        print(f'  3) Number of blades            = {nrblades:6d}')
        print(f'  4) Air density (kg/m^3)        = {rho:10.3f}')
        print(f'  5) Temperature (degr.Celsius)  = {tempc:9.2f}')
        print(f'  6) Prop angle of attack (degr.)= {alphap:9.2f}')
        print(f'  7) Airspeed (m/s)              = {tas:9.2f}')
        print(f'  8) Rpm                         = {rpm:7.0f}')
        print(f'  9) Advance ratio               = {advratio:10.3f}')
        print(f' 10) Scale factor Cd             = {scale_cd:10.3f}')
        print('  ')
        print('  0) Start calculations')
        print('  -1) Quit')

        try:
            # item = int(input('  Select item:'))
            item = 0
        except ValueError:
            print("Invalid input. Please enter an integer.")
            continue

        if item not in range(-1, 11):
            continue

        if item == -1:
            exit()
        elif item == 0:
            break  # Start calculations
        elif item == 1:
            diam = float(input(' Give diameter of prop (m): '))
        elif item == 2:
            beta75 = float(input(' Give beta at 75%: '))
        elif item == 3:
            nrblades = int(input(' Give number of blades: '))
        elif item == 4:
            rho = float(input(' Give air density (kg/m^3): '))
        elif item == 5:
            tempc = float(input(' Give temperature (degr.Celsius): '))
        elif item == 6:
            alphap = float(input(' Give prop angle of attack (degr.): '))
        elif item == 7:
            tas = float(input(' Give airspeed (m/s): '))
            if rpm != 0:
                advratio = tas / ((rpm / 60.0) * diam)
            if advratio != 0:
                rpm = 60.0 * tas / (advratio * diam)
        elif item == 8:
            rpm = float(input(' Give rpm: '))
            if tas != 0:
                advratio = tas / ((rpm / 60.0) * diam)
            if advratio != 0:
                tas = advratio * (rpm / 60.0) * diam
        elif item == 9:
            advratio = float(input(' Give advance ratio: '))
            if tas != 0:
                rpm = 60.0 * tas / (advratio * diam)
            if rpm != 0:
                tas = advratio * (rpm / 60.0) * diam
        elif item == 10:
            scale_cd = float(input(' Enter scale factor for Cd: '))

        if (diam < zero or nrblades == 0 or advratio < zero or
                rho < zero or tas < zero or rpm == 0):
            continue

    # Calculate beta distribution
    beta = beta0 + beta75

    # Calculations
    proparea = (pi / 4) * (diam ** 2)  # m^2
    alphar = alphap * pi / 180.0  # propeller a.o.a. in radians
    tempk = tempc + 273.15  # Kelvin
    q = 0.5 * rho * (tas ** 2)
    omega = 2.0 * pi * (rpm / 60.0)

    # Use Sutherland's law to calculate the viscosity
    mhu = mhuzero * ((tempk / tempkzero) ** 1.5) * (tempkzero + 110) / (tempk + 110)

    # Determine output device (screen or file)
    while True:
        print('\n Select output device    (0) screen')
        print('                         (1) file')
        choice = 1
        try:
            # choice = int(input(' your choice :'))
            choice = 1
        except ValueError:
            print("Invalid input. Please enter 0 or 1.")
            continue

        if choice not in [0, 1]:
            continue

        if choice == 0:
            unit = 6  # Standard output (screen)
            comment = ' '  # no comment when print to screen
            break
        else:
            unit = 222
            while True:
                name2 = os.path.join(output_folder, f"BEM_{prop_name}_output.txt")
                if os.path.exists(name2):
                    # answer = input('"This file exists ! Overwrite (y/n):')
                    answer = 'y'
                    if answer == 'y':
                        try:
                            f_out = open(name2, 'w')  # Open in write mode to overwrite
                            break
                        except Exception as e:
                            print(f"Error opening file: {e}")
                            continue
                    else:
                        continue
                else:
                    try:
                        f_out = open(name2, 'w')  # Open in write mode for new file
                        break
                    except Exception as e:
                        print(f"Error opening file: {e}")
                        continue
            comment = f"Output data {propid}"
            #comment = input(' Enter comment on output file (70 char):')
            break

    # ---------------------------------------------------------------------
    # Now we have all input data, start the calculations
    # ---------------------------------------------------------------------

    # Initialize variables
    dps = 10.0
    psinumber = int(360.0 / dps)  # Number of azimuthal stations
    if psinumber != 12:
        psinumber = 1  # Correct the number if dps is not 30
        dps = 360.0  # if psinumber is 1, then dps is 360 !
    psir = dps * pi / 180.0

    vau0 = np.zeros((12, 20))
    vtu0 = np.zeros((12, 20))
    dptq = np.zeros((12, 20))
    qblade = np.zeros(12)
    tblade = np.zeros(12)
    nblade = np.zeros(12)
    s = np.zeros(12)
    veff_over_v = np.zeros((12, 20))
    gamma = np.zeros((12, 20))

    dptoverq = 0.0
    dndr = 0.0
    dtdr = 0.0
    dqdr = 0.0
    qtotaver = 0.0
    ntotaver = 0.0
    ttotaver = 0.0
    qbladeaver = 0.0
    tbladeaver = 0.0
    nbladeaver = 0.0

    # Print heading
    if choice == 0:
        print(" {:^15} {:^8} {:^8} {:^8} {:^8} {:^8} {:^8} {:^8}".format(
            "r/R", "beta", "phi", "alpha", "Cl", "Cd", "a", "a'"
        ))
    else:
        print(f"{comment:70}", file=f_out)  # Output comment to file
        print(file=f_out)
        print(f" Propeller identification: {propid}", file=f_out)
        print(file=f_out)
        print(
            " {:^15} {:^8} {:^8} {:^8} {:^8} {:^8} {:^8} {:^8}".format(
                "r/R", "beta", "phi", "alpha", "Cl", "Cd", "a", "a'"
            ), file=f_out
        )

    # Loop over blade sections
    for i in range(nrb_sections):
        r = rR[i] * (diam / 2.0)
        dptdrbefore = dptoverq
        dndrbefore = dndr
        dtdrbefore = dtdr
        dqdrbefore = dqdr

        # Loop over azimuth angle psi
        for l in range(psinumber):
            psi = l * dps * pi / 180.0

            # Calculate induced velocity
            if use_induced:
                # Interpolate induced velocities at r/R and psi
                x = psi
                y = rR[i]

                # Find the interval that contains (x, y)
                for m in range(1, nr_psi_i):
                    if psi_i[m - 1] <= x <= psi_i[m]:
                        teta1 = (x - psi_i[m - 1]) / (psi_i[m] - psi_i[m - 1])
                        break
                else:
                    teta1 = 0.0  # Default value if not found

                for m in range(1, nr_rr_i):
                    if rr_i[m - 1] <= y <= rr_i[m]:
                        teta2 = (y - rr_i[m - 1]) / (rr_i[m] - rr_i[m - 1])
                        break
                else:
                    teta2 = 0.0  # Default value if not found

                # Bilinear interpolation of induced velocities
                vx_ii = (1.0 - teta1) * (1.0 - teta2) * vx_i[0, 0] + \
                        teta1 * (1.0 - teta2) * vx_i[1, 0] + \
                        (1.0 - teta1) * teta2 * vx_i[0, 1] + \
                        teta1 * teta2 * vx_i[1, 1]
                vy_ii = (1.0 - teta1) * (1.0 - teta2) * vy_i[0, 0] + \
                        teta1 * (1.0 - teta2) * vy_i[1, 0] + \
                        (1.0 - teta1) * teta2 * vy_i[0, 1] + \
                        teta1 * teta2 * vy_i[1, 1]
                vz_ii = (1.0 - teta1) * (1.0 - teta2) * vz_i[0, 0] + \
                        teta1 * (1.0 - teta2) * vz_i[1, 0] + \
                        (1.0 - teta1) * teta2 * vz_i[0, 1] + \
                        teta1 * teta2 * vz_i[1, 1]

                # Calculate effective velocities
                vplane = tas * np.cos(alphar) + vx_ii
                print(f"vplane: {vplane}")
                vswirl_x = tas * np.sin(alphar) + vy_ii
                vswirl_y = vz_ii
            else:
                vplane = tas * np.cos(alphar)
                print(f"vplane: {vplane}")
                vswirl_x = tas * np.sin(alphar)
                vswirl_y = 0.0

            # Initial estimate for induced velocities
            a_a = 0.0
            a_t = 0.0

            # Iterative calculation of induced velocities
            for d in range(1, 31):
                # Calculate flow angle phi
                print(f"a_a: {a_a}, a_t: {a_t}, omega: {omega}, r: {r}")
                phi = np.arctan(((1.0 + a_a) * vplane) / ((1.0 - a_t) * omega * r))
                print(f"phi: {phi}")
                # Calculate angle of attack alpha
                alpha = np.degrees(phi) - beta[i]
                print(f"alpha: {alpha}")
                # Determine Cl and Cd by interpolation
                re = rho * np.sqrt(vplane ** 2 + (omega * r) ** 2) * cr[i] / mhu
                machi = np.sqrt(vplane ** 2 + (omega * r) ** 2) / 340.0

                # Select correct airfoil data (Re or Mach number)
                if re1[i] > 0 and re2[i] > 0:
                    # Reynolds number is used
                    if re <= (re1[i] + 0.00001):
                        cli = np.interp(alpha, alfa[i, :nr_alfa[i]], cl1[i, :nr_alfa[i]])
                        cdi = np.interp(alpha, alfa[i, :nr_alfa[i]], cd1[i, :nr_alfa[i]])
                    elif re >= (re2[i] - 0.00001):
                        cli = np.interp(alpha, alfa[i, :nr_alfa[i]], cl2[i, :nr_alfa[i]])
                        cdi = np.interp(alpha, alfa[i, :nr_alfa[i]], cd2[i, :nr_alfa[i]])
                    else:
                        # Interpolate between Re1 and Re2
                        factor = (re - re1[i]) / (re2[i] - re1[i])
                        cli = np.interp(alpha, alfa[i, :nr_alfa[i]], cl1[i, :nr_alfa[i]]) + \
                              factor * (np.interp(alpha, alfa[i, :nr_alfa[i]], cl2[i, :nr_alfa[i]]) - \
                                        np.interp(alpha, alfa[i, :nr_alfa[i]], cl1[i, :nr_alfa[i]]))
                        cdi = np.interp(alpha, alfa[i, :nr_alfa[i]], cd1[i, :nr_alfa[i]]) + \
                              factor * (np.interp(alpha, alfa[i, :nr_alfa[i]], cd2[i, :nr_alfa[i]]) - \
                                        np.interp(alpha, alfa[i, :nr_alfa[i]], cd1[i, :nr_alfa[i]]))

                elif mach1[i] > 0 and mach2[i] > 0:
                    # Mach number is used
                    if machi <= (mach1[i] + 0.00001):
                        cli = np.interp(alpha, alfa[i, :nr_alfa[i]], cl1[i, :nr_alfa[i]])
                        cdi = np.interp(alpha, alfa[i, :nr_alfa[i]], cd1[i, :nr_alfa[i]])
                    elif machi >= (mach2[i] - 0.00001):
                        cli = np.interp(alpha, alfa[i, :nr_alfa[i]], cl2[i, :nr_alfa[i]])
                        cdi = np.interp(alpha, alfa[i, :nr_alfa[i]], cd2[i, :nr_alfa[i]])
                    else:
                        # Interpolate between Mach1 and Mach2
                        factor = (machi - mach1[i]) / (mach2[i] - mach1[i])
                        cli = np.interp(alpha, alfa[i, :nr_alfa[i]], cl1[i, :nr_alfa[i]]) + \
                              factor * (np.interp(alpha, alfa[i, :nr_alfa[i]], cl2[i, :nr_alfa[i]]) - \
                                        np.interp(alpha, alfa[i, :nr_alfa[i]], cl1[i, :nr_alfa[i]]))
                        cdi = np.interp(alpha, alfa[i, :nr_alfa[i]], cd1[i, :nr_alfa[i]]) + \
                              factor * (np.interp(alpha, alfa[i, :nr_alfa[i]], cd2[i, :nr_alfa[i]]) - \
                                        np.interp(alpha, alfa[i, :nr_alfa[i]], cd1[i, :nr_alfa[i]]))
                else:
                    # No airfoil data available
                    cli = 0.0
                    cdi = 0.0

                # Scale Cd
                cdi = scale_cd * cdi
                print(f"phi deg: {np.degrees(phi)}")
                # Prandtl's tip loss factor
                f = (nrblades / 2.0) * ((1.0 - rR[i]) / (rR[i] * np.cos(phi)))
                print(f"f: {f}")
                print(f"phi: {phi}")
                if f > 20:
                    fr = 0.0
                else:
                    fr = (2.0 / pi) * np.arccos(np.exp(-f))
                fra = fr

                # Calculate new induced velocities
                sigma = nrblades * cr[i] / (2.0 * pi * r)
                a_a_new = sigma * cli * np.cos(phi) / (4.0 * fra * np.sin(phi) ** 2)
                a_t_new = sigma * cli * np.sin(phi) / (4.0 * fra * np.sin(phi) * np.cos(phi))

                # Relaxation
                a_a = 0.75 * a_a + 0.25 * a_a_new
                a_t = 0.75 * a_t + 0.25 * a_t_new

                # Convergence check
                if abs(a_a - a_a_new) < 0.0001 and abs(a_t - a_t_new) < 0.0001:
                    break

            # Store intermediate results
            vau0[l, i] = a_a
            vtu0[l, i] = a_t
            veff_over_v[l, i] = np.sqrt((1 + a_a) ** 2 * (vplane / tas) ** 2 + (1 - a_t) ** 2 * (
                        omega * r / tas) ** 2)
            gamma[l, i] = 0.5 * cli * cr[i] * veff_over_v[l, i] * tas
            dptq[l, i] = cli * np.cos(phi) * veff_over_v[l, i] ** 2 * (cr[i] / diam) * (nrblades / (2 * pi))
            qblade[l] = qblade[l] + dptq[l, i] * (rR[i] ** 2)
            tblade[l] = tblade[l] + dptq[l, i] * rR[i]
            nblade[l] = nblade[l] + dptq[l, i]

            # Output results for each section
            if choice == 0:
                print(" {:15.4f} {:8.2f} {:8.2f} {:8.2f} {:8.3f} {:8.4f} {:8.3f} {:8.3f}".format(
                    rR[i], beta[i], np.degrees(phi), alpha, cli, cdi, a_a, a_t
                ))
            else:
                print(" {:15.4f} {:8.2f} {:8.2f} {:8.2f} {:8.3f} {:8.4f} {:8.3f} {:8.3f}".format(
                    rR[i], beta[i], np.degrees(phi), alpha, cli, cdi, a_a, a_t
                ), file=f_out)

        # Calculate performance parameters
        dptoverq = dptoverq + 2.0 * pi * rR[i] * dptq[0, i] * (diam / nrb_sections)
        dndr = dndr + 2.0 * pi * rR[i] * dptq[0, i] * tas * (diam / nrb_sections) / omega
        dtdr = dtdr + 2.0 * pi * (rR[i] ** 2) * dptq[0, i] * tas * (diam / nrb_sections)
        dqdr = dqdr + 2.0 * pi * (rR[i] ** 3) * dptq[0, i] * (tas ** 2) * (diam / nrb_sections)

    # End of blade section loop

    # Total performance parameters
    ct = dptoverq
    if ct < zero:
        ct = zero
    cp = omega * dndr / (rho * (tas ** 2))
    if cp < zero:
        cp = zero
    etha = ct / cp * advratio
    power = cp * rho * (tas ** 3) * (diam ** 2)

    # Averaged parameters
    for l in range(psinumber):
        qtotaver = qtotaver + qblade[l] / psinumber
        ntotaver = ntotaver + nblade[l] / psinumber
        ttotaver = ttotaver + tblade[l] / psinumber
        qbladeaver = qbladeaver + qblade[l] / psinumber
        tbladeaver = tbladeaver + tblade[l] / psinumber
        nbladeaver = nbladeaver + nblade[l] / psinumber
        ttot = ttotaver
        qtot = qtotaver

    # Correct torque and thrust for wing incidence
    cnp = ct * np.cos(alphar) + cp / advratio * np.sin(alphar)
    qc = qtot * np.cos(alphar)
    tc = ttot * np.cos(alphar)
    ttot_force = ct * q * proparea ** 2

    qtot_force = cp * q * proparea ** 2 * (diam / 2)

    # Output final results
    if choice == 0:
        print("\n Final results:")
        print(" Propeller identification: ", propid)
        print(" Thrust coefficient        =", ct)
        print(" Power  coefficient        =", cp)
        print(" Efficiency                =", etha)
        print(" Power                     =", power)
        print(f" Total Thrust (ttot): {ttot_force}")
        print(f" Total Torque (qtot): {qtot_force}")
    else:
        print(file=f_out)
        print(" Final results:", file=f_out)
        print(" Propeller identification: ", propid, file=f_out)
        print(" Thrust coefficient        =", ct, file=f_out)
        print(" Power  coefficient        =", cp, file=f_out)
        print(" Efficiency                =", etha, file=f_out)
        print(" Power                     =", power, file=f_out)
        print(" Total Thrust              =", ttot_force, file=f_out)
        print(" Total Torque              =", qtot_force, file=f_out)

    # Close output file
    if choice == 1:
        f_out.close()

    return


propeller_analysis("Hamilton568F", "n\n", "n\n")
