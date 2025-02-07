import os
import subprocess

input_folder = r"C:\Users\tomva\pythonProject\DUUC\analysis_modules"
output_folder = r"C:\Users\tomva\pythonProject\DUUC\data\Polars"
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

    return print(f"\n----- Xfoil output for {file_name}\n"
                 f"Airfoil polar has been created for {airfoil_name}, from "
                 f"{alpha_i} to {alpha_f} in {alpha_step} steps for Reynolds "
                 f"number {Re} with Mach {mach} and saved in {file_name}.txt")

""" Test test if wanted - outdated"""
# xfoil_polar("naca0012", 0, 5, 1, 10000, "test", 0.4)
