import os
import subprocess
from general_functions import print_coordinates, read_file

""" File path definitions"""
input_folder = r"C:\Users\tomva\pythonProject\DUUC\analysis_modules\input_files"
output_folder = r"C:\Users\tomva\pythonProject\DUUC\data\Polars"
xfoil_path = r"C:\Users\tomva\pythonProject\DUUC\analysis_modules\dfdc.exe"
airfoil_coordinates_path = r"C:\Users\tomva\pythonProject\DUUC\data\airfoil_coordinates"
nacelle_coordinates_path = r"C:\Users\tomva\pythonProject\DUUC\data"
dfdc_path = r"C:\Users\tomva\pythonProject\DUUC\analysis_modules\DFDC070\bin\dfdc.exe"

""" definition of the case file"""


def write_case_file(case_name, airfoil_duct, v_inf, v_ref, rpm, rpm2, rho, vso, rmu, altitude, xd_wake,
                    n_wake, lwkrlx, aero_sections, xi_section, a0_deg, dclda, cl_max, cl_min,
                    dclda_stall, dcl_stall, cm_const, m_crit, cd_min, clcd_min, dcddcl, Re_ref,
                    Re_exp, x_disk, n_blades, NRsta, stations, points):
    os.makedirs(input_folder, exist_ok=True)
    case_file_path = os.path.join(input_folder, f"{case_name}.case")

    coordinates_duct = print_coordinates(airfoil_coordinates_path, f"Naca{airfoil_duct}.txt")
    coordinates_nacelle = print_coordinates(nacelle_coordinates_path, f"nacelle_dfdc.txt")
    data_rotor = read_file(nacelle_coordinates_path, "rotor_dfdc.txt")
    drag_data = read_file(nacelle_coordinates_path, "drag_dfdc.txt")

    with open(case_file_path, "w") as file:
        file.write("DFDC Version 0.70\n")
        file.write(f"{case_name}\n\n")
        file.write("OPER\n")
        file.write("!        Vinf         Vref          RPM          RPM\n")
        file.write(f"  {v_inf:<12}   {v_ref:<10}   {rpm:<10}   {rpm2:<10}\n")
        file.write("!         Rho          Vso          Rmu           Alt\n")
        file.write(f"   {rho:<10}   {vso:<10}   {rmu:<12}   {altitude:<12}\n")
        file.write("!       XDwake        Nwake\n")
        file.write(f"  {xd_wake:<10}   {n_wake:<10}\n")
        file.write("!       Lwkrlx\n")
        file.write(f"            {lwkrlx}\n")
        file.write("ENDOPER\n\n")

        file.write("AERO\n")
        file.write("!  #sections\n")
        file.write(f"     {aero_sections}\n")
        file.write("!   Xisection\n")
        file.write(f"    {xi_section}\n")
        file.write("!       A0deg        dCLdA        CLmax         CLmin\n")
        file.write(f"  {a0_deg}         {dclda}        {cl_max}     {cl_min}    \n")
        file.write("!  dCLdAstall     dCLstall      Cmconst         Mcrit\n")
        file.write(f"  {dclda_stall}  {dcl_stall}   {cm_const}      {m_crit}   \n")
        file.write("!       CDmin      CLCDmin     dCDdCL^2\n")
        file.write(f"       {cd_min}    {clcd_min}  {dcddcl}\n")
        file.write("!       REref        REexp\n")
        file.write(f"        {Re_ref}     {Re_exp}    \n")
        file.write("ENDAERO\n\n")

        file.write("ROTOR\n")
        file.write("!       Xdisk        Nblds       NRsta\n")
        file.write(f"       {x_disk}     {n_blades}    {NRsta}\n")
        file.write("!  #stations\n")
        file.write(f"    {stations}\n")
        file.write("!           r        Chord         Beta\n")

        for i in range(len(data_rotor)):
            file.write(f"{data_rotor[i]}\n")

        file.write("ENDROTOR\n\n")

        file.write("DRAGOBJ\n")
        file.write("!  #pts\n")
        file.write(f"   {points}\n")
        file.write("!           x            r          CDA\n")

        for i in range(len(drag_data)):
            file.write(f"{drag_data[i]}\n")

        file.write("ENDDRAGOBJ\n\n")

        file.write("GEOM\n")
        file.write(f"{case_name}\n")

        for i in range(len(coordinates_duct)):
            file.write(f"{coordinates_duct[i]}\n")

        file.write("  999.0 999.0\n")

        for i in range(len(coordinates_nacelle)):
            file.write(f"{coordinates_nacelle[i]}\n")

        file.write("ENDGEOM\n")

        file.close()


def dfdc_operation(file_name, case_name):
    os.makedirs(input_folder, exist_ok=True)
    input_file_path = os.path.join(input_folder, f"{file_name}_input.in")
    case_path = os.path.join(input_folder, f"{case_name}.case")

    with open(input_file_path, "w") as file:
        file.write("LOAD\n")
        file.write(f"{case_path}\n")
        file.write("oper\n")
        file.write("xy\n")
        file.write("exec\n")
        file.write("quit\n")
        file.close()


def run_dfdc(file_name):
    input_file_path = os.path.join(input_folder, f"{file_name}_input.in")

    with open(input_file_path, "r") as f:
        process = subprocess.Popen([dfdc_path], stdin=subprocess.PIPE,
                                   text=True)
        process.communicate(input=f.read())



