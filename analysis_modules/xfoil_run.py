import os
import subprocess
import numpy as np

input_folder = r"C:\Users\tomva\pythonProject\DUUC\analysis_modules\input_files"
output_folder = r"C:\Users\tomva\pythonProject\DUUC\data\Polars"
polar_folder = r"..\data\airfoil_coordinates"
xfoil_path = r"C:\Users\tomva\pythonProject\DUUC\analysis_modules\xfoil.exe"


def xfoil_polar(airfoil_name, alpha_i, alpha_f, alpha_step, Re, file_name, mach):
    os.makedirs(input_folder, exist_ok=True)
    input_file_path = os.path.join(input_folder, f"{file_name}_input.in")
    output_file_path = os.path.join(output_folder, f"{file_name}.txt")

    if os.path.exists(output_file_path):
        os.remove(output_file_path)

    with open(input_file_path, "w") as file:
        file.write(f"{airfoil_name}\n")
        file.write("Oper\n")
        file.write("v\n")
        file.write(f"{Re}\n")
        file.write("mach\n")
        file.write(f"{mach}\n")
        file.write("PACC\n")
        file.write(f"{output_file_path}\n\n")
        file.write("ASeq\n")
        file.write(f"{alpha_i}\n")
        file.write(f"{alpha_f}\n")
        file.write(f"{alpha_step}\n")
        file.write("\n\n")
        file.write("quit\n")
        file.close()

    with open(input_file_path, "r") as f:
        process = subprocess.Popen([xfoil_path], stdin=subprocess.PIPE,
                                   stdout=subprocess.DEVNULL,
                                   text=True)
        process.communicate(input=f.read())
    # old print
    """ print(f"\n----- Xfoil output for {file_name}\n"
          f"Airfoil polar has been created for {airfoil_name}, from "
          f"{alpha_i} to {alpha_f} in {alpha_step} steps for Reynolds "
          f"number {np.round(Re)} with Mach {np.round(mach, 2)} "
          f"and saved in {file_name}.txt")"""
    return print(f"Success: New polar created for: {file_name}")


def xfoil_polar_5seriesload(airfoil_name, alpha_i, alpha_f, alpha_step, Re, file_name, mach):
    os.makedirs(input_folder, exist_ok=True)
    input_file_path = os.path.join(input_folder, f"{file_name}_input.in")
    coordinates_file_path = os.path.join(polar_folder, f"{airfoil_name}")
    output_file_path = os.path.join(output_folder, f"{file_name}.txt")

    if os.path.exists(output_file_path):
        os.remove(output_file_path)

    with open(input_file_path, "w") as file:
        file.write(f"load\n")
        file.write(f"{coordinates_file_path}.txt\n")
        file.write(f"naca43013\n")
        file.write("Oper\n")
        file.write("v\n")
        file.write(f"{Re}\n")
        file.write("mach\n")
        file.write(f"{mach}\n")
        file.write("PACC\n")
        file.write(f"{output_file_path}\n\n")
        file.write("ASeq\n")
        file.write(f"{alpha_i}\n")
        file.write(f"{alpha_f}\n")
        file.write(f"{alpha_step}\n")
        file.write("\n\n")
        file.write("quit\n")
        file.close()

    with open(input_file_path, "r") as f:
        process = subprocess.Popen([xfoil_path], stdin=subprocess.DEVNULL,
                                   stdout=subprocess.DEVNULL,
                                   text=True)
        process.communicate(input=f.read())

    return print(f"\n----- Xfoil output for {file_name}\n"
                 f"Airfoil polar has been created for {airfoil_name}, from "
                 f"{alpha_i} to {alpha_f} in {alpha_step} steps for Reynolds "
                 f"number {np.round(Re)} with Mach {np.round(mach, 2)} "
                 f"and saved in {file_name}.txt")


def create_new_polars(profile_pylon, profile_support, profile_duct, profile_control, re_pylon, re_support, re_duct,
                      re_control, mach):
    print("XFoil starts......")
    airfoil_pylon = "Naca" + profile_pylon
    file_pylon = "pylon" + profile_pylon
    airfoil_support = "Naca" + profile_support
    file_support = "support" + profile_support
    airfoil_duct = "Naca" + profile_duct
    file_duct = "duct" + profile_duct
    airfoil_control = "Naca" + profile_control
    file_control = "control" + profile_control

    airfoils = [airfoil_pylon, airfoil_support, airfoil_duct, airfoil_control]
    filenames = [file_pylon, file_support, file_duct, file_control]
    re_comp = [re_pylon, re_support, re_duct, re_control]

    for i in range(len(airfoils)):
        xfoil_polar(airfoils[i], -5, 15, 1, re_comp[i], filenames[i], mach)
    print("XFoil closes......")

""" Test test if wanted"""

xfoil_polar("Naca0016", -6, 20, 0.5, 4.698e6, "support0016", 0.44)
