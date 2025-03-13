import os
import numpy as np
import matplotlib.pyplot as plt


def print_coordinates(file_path, file_name: str):
    matrix = []
    file_path = os.path.join(file_path, f"{file_name}")
    with open(file_path, "r") as coordinates:
        for line in coordinates:
            matrix.append(line.strip())
    return matrix


def write_coordinates_nacelle(nacelle_length, nacelle_diameter, file_name="nacelle_dfdc.txt"):
    """ Assume that the nacelle is an axi_symmetric body with a front and end rounding"""
    write_to_filepath = r"C:\Users\tomva\pythonProject\DUUC\data"
    file_path = os.path.join(write_to_filepath, file_name)

    radius = 0.5 * nacelle_diameter

    coordinates = []

    # First quarter circle (front part)
    theta = np.linspace(np.pi, np.pi / 2, 20)  # 20 points for smooth curve
    for angle in theta:
        x = radius * np.cos(angle) + radius
        y = radius * np.sin(angle)
        coordinates.append((x, y))

    # Straight section
    x_values = np.linspace(radius, nacelle_length - radius, 30)  # 30 points for the straight part
    for x in x_values:
        coordinates.append((x, radius))

    # Final quarter circle (back part)
    theta = np.linspace(np.pi / 2, 0, 20)  # 20 points for smooth curve
    for angle in theta:
        x = (nacelle_length - radius) + radius * np.cos(angle)
        y = radius * np.sin(angle)
        coordinates.append((x, y))

    # Write to file
    with open(file_path, "w") as file:
        for x, y in coordinates:
            file.write(f"{x:.6f} {y:.6f}\n")


def read_file(file_path, file_name: str):
    matrix = []
    file_path = os.path.join(file_path, f"{file_name}")
    with open(file_path, "r") as lines:
        for line in lines:
            matrix.append(line.strip())
    return matrix


def shp_conversion(shp, unit: str):
    if unit == "SI":
        power = shp * 745.7  # output in [W]
        return power
    if unit == "metric":
        power = shp * 0.986  # output in [hp]
        return power
    else:
        print("Wrong unit input. Options: SI or metric")
        return None


def write_blade_data(prop_name, prop_diam, n_blades, temperature, density, alpha, tas, advance,
                     beta75, n_sections):
    input_folder = r"C:\Users\tomva\pythonProject\DUUC\data"
    geometry_folder = r"C:\Users\tomva\pythonProject\DUUC\data\airfoil_coordinates"
    os.makedirs(input_folder, exist_ok=True)
    input_file_path = os.path.join(input_folder, f"{prop_name}_data.txt")

    with open(input_file_path, "w") as file:
        file.write(f"{prop_name}\n\n")
        file.write(f"{prop_diam}\n")
        file.write(f"{n_blades}\n")
        file.write(f"{temperature}\n")
        file.write(f"{density}\n")
        file.write(f"{alpha}\n")
        file.write(f"{tas}\n")
        file.write(f"{advance}\n")
        file.write(f"{beta75}\n")
        file.write(f"{n_sections}\n")

        file_path = os.path.join(geometry_folder, f"{prop_name}_geometry.txt")
        with open(file_path, "r") as coordinates:
            next(coordinates)
            for line in coordinates:
                file.write(f"{line}")

        file.write("same\n")
        file.write(f"{prop_name}_airfoil.txt")
        file.close()
    return print(f"file written")
