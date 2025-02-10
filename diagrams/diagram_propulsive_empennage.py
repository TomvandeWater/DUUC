import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# load airfoil coordinates
file_path = r'/data/Naca0012.txt'

start_line1 = 0
data_cp_init = np.loadtxt(file_path)
x_airfoil = data_cp_init[:, 0]
y_airfoil = data_cp_init[:, 1]

# Create a figure and an axis
fig, ax = plt.subplots()

# Pylon
ax.vlines(0.6, 0, 0.3-0.058, colors='black', linewidth=1.5)
ax.vlines(1.0, 0, 0.3-0.03, colors='black', linewidth=1.5)

# lower airfoil
plt.plot(x_airfoil + 0.2, y_airfoil + 0.3, color='black')

# upper airfoil
plt.plot(x_airfoil + 0.2, y_airfoil + 2.2, color='black')

# nacelle
nacelle = patches.Rectangle((0.5, 1.15), 0.7, 0.2, linewidth=1.5,
                            edgecolor='black', facecolor='none')
ax.add_patch(nacelle)
half_circle_nacelle = patches.Arc((1.2, 1.25), 2*0.1, 2*0.1, angle=-90,
                                  theta1=0, theta2=180, color='black',
                                  linewidth=1.5)
ax.add_patch(half_circle_nacelle)

# propeller
rounded_square = patches.FancyBboxPatch((0.46, 0.43), 0.005, 1.65,
                                        boxstyle="round,pad=0.03,rounding_size="+str(0.05),
                                        linewidth=1.5, edgecolor='black',
                                        facecolor='none')
half_circle_prop = patches.Arc((0.43, 1.25), 2*0.1, 2*0.1, angle=90, theta1=0,
                               theta2=180, color='black', linewidth=1.5)
ax.add_patch(half_circle_prop)

# support
ax.vlines(0.6, 1.35, 2.2-0.058, colors='black', linewidth=1.5)
ax.vlines(0.9, 1.35, 2.2-0.03, colors='black', linewidth=1.5)

ax.vlines(0.6, 0.3+0.058, 1.15, colors='black', linewidth=1.5)
ax.vlines(0.9, 0.3+0.03, 1.15, colors='black', linewidth=1.5)

# control vanes
hcv1 = patches.Rectangle((1.10, 0.33), 0.1, 0.82, linewidth=1.5,
                         edgecolor='black', facecolor='none')
ax.add_patch(hcv1)
hcv2 = patches.Rectangle((1.10, 1.35), 0.1, 0.82, linewidth=1.5,
                         edgecolor='black', facecolor='none')
ax.add_patch(hcv2)

plt.plot(x_airfoil*0.1 + 1.1, y_airfoil*0.3 + 1.25, color='black')


# Add the square to the plot
ax.add_patch(rounded_square)


# Set the limits of the plot to make sure the square is visible
ax.set_xlim(0, 1.5)
ax.set_ylim(0, 2.5)

# Display the plot
plt.show()
