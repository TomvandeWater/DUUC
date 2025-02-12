import os
import numpy as np
import subprocess
import constants
import data.atr_reference as ref
import flow_conditions


""" file paths"""
input_folder = r"C:\Users\tomva\pythonProject\DUUC\analysis_modules\input_files"
output_folder = r"C:\Users\tomva\pythonProject\DUUC\data\Polars"
avl_path = r"C:\Users\tomva\pythonProject\DUUC\analysis_modules\avl.exe"


def write_geometry(file_name, mach, profile_vt: str, profile_ht: str, s_ref, c_ref, b_ref, x_ref,
                   y_ref, z_ref, cd0_ref, c_root_ht, tr_ht, sweep_ht, b_ht,
                   x_le_vt, c_root_vt, tr_vt, sweep_vt, b_vt):
    """ Note that this equation is based on a t-tail configuration """
    input_file_path = os.path.join(input_folder, f"{file_name}_input.avl")

    c_vt_tip = c_root_vt * tr_vt
    x_le_vt_tip = b_vt / 2 * np.sin(np.radians(sweep_vt))

    c_ht_tip = c_root_ht * tr_ht
    x_le_ht_tip = b_ht / 2 * np.sin(np.radians(sweep_ht)) + x_le_vt_tip

    with open(input_file_path, "w") as file:
        file.write(f"{file_name}\n")
        file.write(f"{mach}\n")
        file.write(f"0          0           0\n")
        file.write(f"{s_ref}    {c_ref}     {b_ref}\n")
        file.write(f"{x_ref}    {y_ref}     {z_ref}\n")
        file.write(f"{cd0_ref}\n")
        file.write("#\n#===========================================\n#\n")

        file.write(f"SURFACE\n")
        file.write("HORIZONTAL TAIL\n")
        file.write("#Nchordwise  Cspace   Nspanwise   Sspace\n")
        file.write(f"{10}         {1}      {20}        {1}\n")
        file.write("#\nYDUPLICATE\n0.0\n#\nANGLE\n0.00\n#\n# x,y,z bias for whole surface\nTRANSLATE\n")
        file.write("0.00    0.00      0.00\n")
        file.write("#----------------------------------------------------------\n")
        file.write("    Xle         Yle         Zle         chord       angle   Nspan  Sspace\n")
        file.write("SECTION\n")
        file.write(f"{x_le_vt_tip}      {0}       {b_vt}     {c_root_ht}     0.0\n")
        file.write(f"NACA\n{profile_ht}\n\n")
        file.write("#----------------------------------------------------------\n")
        file.write("SECTION\n")
        file.write(f"{x_le_ht_tip}      {b_ht / 2}       {b_vt}     {c_ht_tip}     0.0\n")
        file.write(f"NACA\n{profile_ht}\n\n")

        file.write("#\n#===========================================\n#\n")
        file.write(f"SURFACE\n")
        file.write("VERTICAL TAIL\n")
        file.write("#Nchordwise  Cspace   Nspanwise   Sspace\n")
        file.write(f"{10}         {1}      {20}        {1}\n")
        # file.write("#\nYDUPLICATE\n#\nANGLE\n0.00\n#\n# x,y,z bias for whole surface\nTRANSLATE\n")
        # file.write("0.00    0.00      0.00\n")
        file.write("#----------------------------------------------------------\n")
        file.write("    Xle         Yle         Zle         chord       angle   Nspan  Sspace\n")
        file.write("SECTION\n")
        file.write(f"{x_le_vt}      {0}       {0}     {c_root_vt}     0.0\n")
        file.write(f"NACA\n{profile_vt}\n\n")
        file.write("#----------------------------------------------------------\n")
        file.write("SECTION\n")
        file.write(f"{x_le_vt_tip}      {0}       {b_vt}     {c_vt_tip}     0.0\n")
        file.write(f"NACA\n{profile_vt}\n")
        file.close()
    return print("-- AVL -- Geometry file successfully printed")


def write_mass(file_name, m_vt, m_ht, b_vt, b_ht, c_vt):
    """ Intertias are optional and will not be included for now, point masses are assumed to be
    at half of the surfaces at half of the chord"""
    input_file_path = os.path.join(input_folder, f"{file_name}.mass")

    with open(input_file_path, "w") as file:
        file.write(f"{file_name}\n")
        file.write("#x, y, z coordinate system matches AVL default\n")
        file.write(f"Lunit = 1 m\n")
        file.write(f"Munit = 1 kg\n")
        file.write(f"Tunit = 1 s\n")
        file.write("#----------------------------------------------------------\n")
        file.write(f"g = {constants.g}\n")
        file.write(f"rho = {flow_conditions.rho}\n")
        file.write("#----------------------------------------------------------\n")
        file.write("mass   x   y   z    Ixx    Iyy   Izz    Ixy   Ixz  Iyz\n")
        # vertical tail
        file.write(f"{m_vt}   {c_vt / 2} {0} {b_vt / 2} !Vertical tail\n")
        # horizontal tail right
        file.write(f"{m_ht / 2}   {c_vt / 2} {b_ht / 2} {b_vt} !Horizontal tail\n")
        # horizontal tail left
        file.write(f"{m_ht / 2}   {c_vt / 2} {-b_ht / 2} {b_vt} !Horizontal tail\n")
        file.close()
    return print("-- AVL -- Mass file successfully printed")


def write_case(file_name, alpha, velocity):
    input_file_path = os.path.join(input_folder, f"{file_name}.run")

    with open(input_file_path, "w") as file:
        file.write(f"-------------------------------\n")
        file.write("Run cas  1: case number 1\n")
        # constraints block
        file.write("alpha   -> CL = 0.00\n")
        file.write("beta   -> beta = 0.00\n")
        file.write("pb/2V   -> pb/2V = 0.00\n")
        file.write("qc/2V   -> qc/2V = 0.00\n")
        file.write("rb/2V   -> rb/2V = 0.00\n")
        file.write("camber   -> camber = 0.00\n")
        file.write("aileron   -> Cl roll mom = 0.00\n")
        file.write("elevator   -> Cm pitchmom = 0.00\n")
        file.write("rudder   -> Cn yaw mom = 0.0\n\n")
        # flight conditions
        file.write(f"alpha = {alpha}  deg\n")
        file.write(f"beta = 0.00  deg\n")
        file.write(f"pb/2V = 0.00  \n")
        file.write(f"qc/2V = 0.00  \n")
        file.write(f"rb/2V = 0.00  \n")
        file.write(f"CL = 0.00  \n")
        file.write(f"Cdo = 0.00  \n")
        file.write(f"bank = 0.00  deg\n")
        file.write(f"elevation = 0.00  deg\n")
        file.write(f"heading = 0.00  deg\n")
        file.write(f"Mach = {flow_conditions.Mach}\n")
        file.write(f"velocity = {velocity}  m/s\n")
        file.write(f"density = {flow_conditions.rho} kg/m^3\n")
        file.write(f"grav.acc. = {constants.g}  m/s^2\n")
        file.write(f"turn_rad. = 0.00  m\n")
        file.write(f"load_fac. = 0.00  \n")
        # viscous properties
        file.write("visc CL_a = 0.00\n")
        file.write("visc CL_u = 0.00\n")
        file.write("visc CM_a = 0.00\n")
        file.write("visc CM_u = 0.00\n")
        file.close()
    return print("-- AVL -- Case file successfully printed")


def run_avl(file_name, alpha, velocity, mach, ref_matrix):
    input_file = os.path.join(input_folder, f"{file_name}.in")

    # READ DATA FROM REFERENCE DATA MATRIX
    airfoil_ht = ref_matrix[0]
    airfoil_vt = ref_matrix[1]
    c_root_h = ref_matrix[2]
    c_root_v = ref_matrix[3]
    tr_h = ref_matrix[4]
    tr_v = ref_matrix[5]
    phi_qc_h = ref_matrix[6]
    phi_qc_v = ref_matrix[7]
    b_h = ref_matrix[8]
    b_v = ref_matrix[9]
    m_ht = ref_matrix[10]
    m_vt = ref_matrix[11]
    cd0_ref = ref_matrix[12]

    # write geometry, mass and case file
    write_geometry(file_name, mach, airfoil_ht, airfoil_vt, 10, 2, 25, 0,
                   0, 0, cd0_ref, c_root_h, tr_h, phi_qc_h, b_h, 0, c_root_v,
                   tr_v, phi_qc_v, b_v)
    write_mass(file_name, m_vt, m_ht, b_v, b_h, c_root_v)
    write_case(file_name, alpha, velocity)

    with open(input_file, "w") as file:
        file.write(f"LOAD\n")
        file.write(f"C:\Users\tomva\pythonProject\DUUC\analysis_modules\input_files\{file_name}_input.avl\n")


        file.close()

    # run AVL
    with open(input_file, "r") as f:
        process = subprocess.Popen([avl_path], stdin=subprocess.PIPE,
                                   text=True)
        process.communicate(input=f.read())
    return print("-- AVL -- run completed -> output file written")

""" test section"""
"""# write_geometry("test_avl", flow_conditions.Mach, ref.airfoil_ht, ref.airfoil_vt, 10,
               # 2, 25, 0, 0, 0, 0, ref.c_root_h, ref.tr_h,
               # ref.phi_qc_h, ref.b_h, 0, ref.c_root_v, ref.tr_v, ref.phi_qc_v, ref.b_v)
# write_mass("test_avl", ref.m_vt, ref.m_ht, ref.b_v, ref.b_h, ref.c_root_v)
# write_case("test_avl", 10, 181) """


