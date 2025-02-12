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

