import numpy as np
from scipy.interpolate import interp1d
import os
import re

data_folder = r"C:\Users\tomva\pythonProject\DUUC\data"

# -----------------------------------------------------------------------
# Program to calculate propeller flow based on the blade element
# theory.
# Prandtl's tip loss factor is used
#
# Blade elements characteristics are read from file
#
#     first blade section at r/R <= 0.75
#     const. nr. of mach- or Reynolds-numbers    :2
#     max. nr. of angles of attack               :28
#     max. nr. of blade sections                 :10
#     (max. nr. of delta-psi                     :360/30=12 )
#
#     gamma     = non-dimensional local circulation
#
# Author       : L.L.M. Veldhuis
#                Delft University of Technology
#                Faculty of Aerospace Engineering
# Init Date    : november 1989
#
# Adaptations:
#
#   6 may 2001 : psi_number is set to either 1 or 12 gain using if
#                statement
#  15 may 2001 : option to change Cd of blade airfoil with scalefactor.
#   4 oct 2001 : number of blade sections increased
# -----------------------------------------------------------------------


def propeller_analysis():
    """
    Calculates propeller flow based on blade element theory.
    """

    # Constants and default values
    ZERO = 1.0e-20
    PI = 4.0 * np.arctan(1.0)
    TEMP_K_ZERO = 288.16  # at zero S.A
    MHU_ZERO = 1.7894e-5  # idem
    GNU_OUTPUT = 1  # default value for gnu output
    MAX_PSI_I = 30
    MAX_RR_I = 25

    new_case = 0
    use_induced = False
    ducted = False
    scale_cd = 1.0
    gnuplotmacros = '~/surfdrive/5A-UTILITIES-1/fortran/prop1b/gnuplotmacros/'

    # Clear screen
    # os.system('clear')  # UNIX

    if new_case == 0:
        print(" *---------------------------------------------*")
        print(" |                 P R O P 1 B                 |")
        print(" |                                             |")
        print(" |                 version 2.2                 |")
        print(" |                                             |")
        print(" |            Author: L.L.M. Veldhuis          |")
        print(" |                                             |")
        print(" |        Delft University of Technology       |")
        print(" |       Faculty of Aerospace Engineering      |")
        print(" |                                             |")
        print(" *---------------------------------------------*")

    # Get blade data file
    while True:
        blade_data = os.path.join(data_folder, f"blade_data.txt")
        try:
            with open(blade_data, 'r') as f:
                # Read header information
                propid = f.readline().strip()
                f.readline()  # Empty line
                # Handle integer/float conversions safely
                diam = float(f.readline().strip())
                nrblades = int(float(f.readline().strip()))  # Convert to float first, then int
                tempc = float(f.readline().strip())
                rho = float(f.readline().strip())
                alphap = float(f.readline().strip())
                tas = float(f.readline().strip())
                advratio = float(f.readline().strip())
                beta75 = float(f.readline().strip())
                nrb_sections = int(float(f.readline().strip()))  # Same fix here

                # Read blade geometry
                rR = np.zeros(nrb_sections)
                cr = np.zeros(nrb_sections)
                beta0 = np.zeros(nrb_sections)
                for i in range(nrb_sections):
                    line = f.readline().strip()

                    # Handle non-standard whitespace (tabs, multiple spaces)
                    clean_line = line.split('#')[0]
                    parts = re.split(r'\s+', clean_line.strip())  # Split on any whitespace

                    if len(parts) != 3:
                        raise ValueError(f"Line {i + 1}: Need exactly 3 values. Found: '{line}'")

                    rR[i], cr[i], beta0[i] = map(float, parts)

                # Tip r over R should be slightly smaller than 1.0
                if abs(rR[nrb_sections - 1] - 1.0) < ZERO:
                    rR[nrb_sections - 1] = 0.99

                # Read airfoil characteristics
                samediff = f.readline().strip()
                if samediff == 'same':
                    airfoilfile_name = f.readline().strip()
                    airfoilfile = os.path.join(data_folder, f"{airfoilfile_name}")

                nr_alfa = np.zeros(nrb_sections, dtype=int)
                ma1_or_re1 = np.zeros(nrb_sections)
                ma2_or_re2 = np.zeros(nrb_sections)
                alfa = np.zeros((nrb_sections, 28))
                cl1 = np.zeros((nrb_sections, 28))
                cd1 = np.zeros((nrb_sections, 28))
                cl2 = np.zeros((nrb_sections, 28))
                cd2 = np.zeros((nrb_sections, 28))

                for i in range(nrb_sections):
                    if samediff == 'same':
                        af_file = airfoilfile
                    else:
                        af_file = f.readline().strip()
                    with open(af_file, 'r') as af:
                        nr_alfa[i] = int(af.readline())
                        ma1_or_re1[i] = float(af.readline())

                        for j in range(nr_alfa[i]):
                            alfa[i, j], cl1[i, j], cd1[i, j] = map(float, af.readline().split())

                        ma2_or_re2[i] = float(af.readline())
                        for j in range(nr_alfa[i]):
                            alfa[i, j], cl2[i, j], cd2[i, j] = map(float, af.readline().split())

                break  # Exit loop if file reading was successful

        except FileNotFoundError:
            print(f"Error: File not found.")
        except Exception as e:
            print(f"Error reading file: {e}")

    # Calculate RPM based on input
    rpm = (tas / (advratio * diam)) * 60  # RPM is calculated

    # Derived parameters
    beta = beta0 + beta75
    proparea = (PI / 4) * (diam ** 2)  # m^2
    alphar = alphap * PI / 180.0  # propeller a.o.a. in radians
    tempk = tempc + 273.15  # Kelvin
    q = 0.5 * rho * (tas ** 2)
    omega = 2.0 * PI * (rpm / 60.0)

    # Sutherland's law for viscosity
    mhu = MHU_ZERO * ((tempk / TEMP_K_ZERO) ** 1.5) * (TEMP_K_ZERO + 110) / (tempk + 110)

    # Ducted propeller?
    answer = input("Ducted propeller? (y/n): ")
    ducted = answer.lower() == 'y'

    # Use wing induced velocity data?
    answer = input("Use (wing) induced velocity data? (y/n): ")
    use_induced = answer.lower() == 'y'

    if use_induced:
        file_induced = input("Enter file containing induced velocities: ")

        # Check file contents
        while True:
            try:
                ifilecontents = int(input("Select file contents:\n"
                                          " (1) psi,r/R,Vx,Vy,Vz\n"
                                          " (2) psi,r/R,Vx (Vy=Vz=0)\n"
                                          "Your choice: "))
                if ifilecontents in [1, 2]:
                    break
                else:
                    print("Invalid choice. Please enter 1 or 2.")
            except ValueError:
                print("Invalid input. Please enter a number.")

        try:
            with open(file_induced, 'r') as f:
                nr_psi_i, nr_rr_i = map(int, f.readline().split())

                if nr_psi_i > MAX_PSI_I or nr_rr_i > MAX_RR_I:
                    print("Too many induced data points")
                    return

                psi_i = np.zeros(nr_psi_i)
                rr_i = np.zeros(nr_rr_i)
                vx_i = np.zeros((nr_psi_i, nr_rr_i))
                vy_i = np.zeros((nr_psi_i, nr_rr_i))
                vz_i = np.zeros((nr_psi_i, nr_rr_i))

                if ifilecontents == 1:
                    for i in range(nr_psi_i):
                        for j in range(nr_rr_i):
                            psi_i[i], rr_i[j], vx_i[i, j], vy_i[i, j], vz_i[i, j] = map(float, f.readline().split())
                else:
                    for i in range(nr_psi_i):
                        for j in range(nr_rr_i):
                            psi_i[i], rr_i[j], vx_i[i, j] = map(float, f.readline().split())
                            vy_i[i, j] = 0.0
                            vz_i[i, j] = 0.0

        except FileNotFoundError:
            print("Error: Induced velocity file not found.")
            return
        except Exception as e:
            print(f"Error reading induced velocity file: {e}")
            return

    # Interactive input data loop
    while True:
        # os.system('clear')  # UNIX

        print(" Input data:")
        print(f" 1) Diameter of prop (m)        = {diam:10.3f}")
        print(f" 2) Beta at 75% (degr.)         = {beta75:9.2f}")
        print(f" 3) Number of blades            = {nrblades:6d}")
        print(f" 4) Air density (kg/m^3)        = {rho:10.3f}")
        print(f" 5) Temperature (degr.Celsius)  = {tempc:9.2f}")
        print(f" 6) Prop angle of attack (degr.)= {alphap:9.2f}")
        print(f" 7) Airspeed (m/s)              = {tas:9.2f}")
        print(f" 8) Rpm                         = {rpm:7.0f}")
        print(f" 9) Advance ratio               = {advratio:10.3f}")
        print(f"10) Scale factor Cd             = {scale_cd:10.3f}")
        print("")
        print(" 0) Start calculations")
        print("")
        print(" -1) Quit")
        print("")

        while True:
            try:
                item = int(input(" Select item: "))
                break
            except ValueError:
                print("Invalid input. Please enter an integer.")

        if item == -1:
            return

        if item == 0:
            break

        if item == 1:
            diam = float(input("Give diameter of prop (m): "))
        elif item == 2:
            beta75 = float(input("Give beta at 75%: "))
        elif item == 3:
            nrblades = int(input("Give number of blades: "))
        elif item == 4:
            rho = float(input("Give air density (kg/m^3): "))
        elif item == 5:
            tempc = float(input("Give temperature (degr.Celsius): "))
        elif item == 6:
            alphap = float(input("Give prop angle of attack (degr.): "))
        elif item == 7:
            tas = float(input("Give airspeed (m/s): "))
            if rpm != 0:
                advratio = tas / ((rpm / 60.0) * diam)
            if advratio != 0:
                rpm = 60.0 * tas / (advratio * diam)
        elif item == 8:
            rpm = float(input("Give rpm: "))
            if tas != 0:
                advratio = tas / ((rpm / 60.0) * diam)
            if advratio != 0:
                tas = advratio * (rpm / 60.0) * diam
        elif item == 9:
            advratio = float(input("Give advance ratio: "))
            if tas != 0:
                rpm = 60.0 * tas / (advratio * diam)
            if rpm != 0:
                tas = advratio * (rpm / 60.0) * diam
        elif item == 10:
            scale_cd = float(input("Enter scale factor for Cd: "))
        else:
            print("Invalid item selected.")
            continue  # restart the loop

        if (diam < ZERO or nrblades == 0 or advratio < ZERO or
            rho < ZERO or tas < ZERO or rpm == 0):
            continue  # restart the loop

        # Recalculate beta distribution
        beta = beta0 + beta75

        # Recalculate derived parameters
        proparea = (PI / 4) * (diam ** 2)
        alphar = alphap * PI / 180.0
        tempk = tempc + 273.15
        q = 0.5 * rho * (tas ** 2)
        omega = 2.0 * PI * (rpm / 60.0)
        mhu = MHU_ZERO * ((tempk / TEMP_K_ZERO) ** 1.5) * (TEMP_K_ZERO + 110) / (tempk + 110)

    # Determine output device (screen or file)
    while True:
        try:
            choice = int(input("\nSelect output device\n"
                               " (0) screen\n"
                               " (1) file\n"
                               "your choice: "))
            if choice in [0, 1]:
                break
            else:
                print("Invalid choice. Please enter 0 or 1.")
        except ValueError:
            print("Invalid input. Please enter a number.")

    if choice == 0:
        unit = 6
        comment = ' '  # no comment when print to screen
    else:
        unit = 222
        while True:
            name2 = input("Enter output file name: ")
            if os.path.exists(name2):
                answer = input("\"This file exists! Overwrite (y/n): ")
                if answer.lower() == 'y':
                    try:
                        f = open(name2, 'w')  # open in write mode to overwrite
                        break  # if successful, exit the loop
                    except Exception as e:
                        print(f"Error opening file: {e}")
                else:
                    continue  # if not overwriting, ask for the filename again
            else:
                try:
                    f = open(name2, 'w')  # open in write mode to create the file
                    break  # if successful, exit the loop
                except Exception as e:
                    print(f"Error opening file: {e}")

    #-----------------------------------------------------------------------
    # Blade element loop
    #-----------------------------------------------------------------------
    dr=0.001
    nr_segments=int(1/dr)
    dqdr_list=[]
    dtdr_list=[]
    dndr_list=[]
    rR_segments=np.zeros(nr_segments)
    r_segments=np.zeros(nr_segments)

    for i in range(nr_segments):
        rR_segments[i]=(i+1)*dr
        r_segments[i]=rR_segments[i]*diam/2

    a_a=0
    a_t=0
    vplane=tas*np.cos(alphar)
    van=tas*np.sin(alphar)
    cp=0
    ct=0

    for i in range(1,nr_segments):
            r=r_segments[i]
            rR_segment=rR_segments[i]
            # Interpolate blade characteristics

            f_cr = interp1d(rR, cr, kind='linear', fill_value="extrapolate")
            cr_segment = f_cr(rR_segment)

            f_beta = interp1d(rR, beta, kind='linear', fill_value="extrapolate")
            beta_segment = f_beta(rR_segment)

            vt=omega*r
            phi=np.arctan((vplane*(1+a_a))/(vt*(1-a_t)))
            alpha=np.degrees(phi)-beta_segment

            # Interpolate cl and cd
            f_cl1 = interp1d(alfa[0,:nr_alfa[0]], cl1[0,:nr_alfa[0]], kind='linear', fill_value="extrapolate")
            cl1_segment = f_cl1(alpha)

            f_cd1 = interp1d(alfa[0,:nr_alfa[0]], cd1[0,:nr_alfa[0]], kind='linear', fill_value="extrapolate")
            cd1_segment = f_cd1(alpha)*scale_cd

            cli=cl1_segment
            cdi=cd1_segment

            dqdr=(0.5*rho*((vplane*(1+a_a))**2+(vt*(1-a_t))**2)*cr_segment*cli*np.sin(phi))*r*nrblades
            dtdr=(0.5*rho*((vplane*(1+a_a))**2+(vt*(1-a_t))**2)*cr_segment*cli*np.cos(phi))*nrblades
            dndr=(0.5*rho*((vplane*(1+a_a))**2+(vt*(1-a_t))**2)*cr_segment*cdi)*nrblades
            dqdr_list.append(dqdr)
            dtdr_list.append(dtdr)
            dndr_list.append(dndr)

            cp=cp+dqdr*dr
            ct=ct+dtdr*dr

    cp=cp/(rho*(rpm/60)**3*diam**5)
    ct=ct/(rho*(rpm/60)**2*diam**4)

    #Performance factors
    etha=advratio*ct/cp

    print(f"Advance Ratio: {advratio}")
    print(f"Power Coefficient: {cp}")
    print(f"Thrust Coefficient: {ct}")
    print(f"Efficiency: {etha}")



# Run the analysis
propeller_analysis()
