import os
import subprocess


def xfoil_polar(airfoil_name, alpha_i, alpha_f, alpha_step, Re, file_name, mach):
    if os.path.exists("{0}.txt".format(file_name)):
        os.remove("{0}.txt".format(file_name))

    file = open("{0}_input.in".format(file_name), "w")
    file.write("{0}\n".format(airfoil_name))
    file.write("Oper\n")
    file.write("v\n")
    file.write("{0}\n".format(Re))
    file.write("mach\n")
    file.write("{0}\n".format(mach))
    file.write("PACC\n")
    file.write("{0}.txt\n\n".format(file_name))
    file.write("ASeq\n")
    file.write("{0}\n".format(alpha_i))
    file.write("{0}\n".format(alpha_f))
    file.write("{0}\n".format(alpha_step))
    file.write("\n\n")
    file.write("quit\n")
    file.close()

    # subprocess.call("xfoil.exe < input_file.in>", shell=True)
    process = subprocess.Popen(["xfoil.exe"], stdin=subprocess.PIPE,
                               text=True)
    process.communicate(input=open("{0}_input.in".format(file_name)).read())
    return print("Airfoil polar has been created for {0}, from {1} to {2} "
                 "in {3} steps for Reynolds numer {4} and saved in {5}".format(
                  airfoil_name, alpha_i, alpha_f, alpha_step, Re, file_name))


# xfoil_polar("naca0012", 0, 5, 1, 10000, "test", 0.4)
