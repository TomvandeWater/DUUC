import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker


def plot_inflow_properties(v_input, a_input, station):

    v_input2 = [100, 150, 175, 150, 125, 150, 180, 100]
    a_input2 = [0, 10, 20, 30, 25, 10, 25, 0]

    fig, ax1 = plt.subplots(figsize=(12, 8))
    fig.canvas.manager.set_window_title("Inflow parameters per component")
    img = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\propulsive_empennage.png")
    ax1.imshow(img, extent=[0.35, 6.00, 0.00, 200.00])
    ax1.plot(station, v_input, label=r"$V_{inflow}$ - power on", color="tab:blue", marker='o')
    ax1.plot(station, v_input2, label=r"$V_{inflow}$ - power off", color="tab:blue", linestyle="dotted", marker='o')
    ax1.set_xlim(0, 6)
    ax1.set_ylim(0, 200)
    ax1.xaxis.set_major_locator(ticker.NullLocator())
    ax1.set_aspect(aspect=0.02)
    ax1.set_ylabel(r'$v_{\infty}$ [m/s]', color="tab:blue")
    ax1.set_title(r'Inflow parameters per component')
    ax1.tick_params(axis='y', labelcolor="tab:blue")

    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(station, a_input, label=r"$\alpha$ - power on", color="tab:orange", marker='o')
    ax2.plot(station, a_input2, label=r"$\alpha$ - power off", color="tab:orange", linestyle="dotted", marker='o')

    ax2.set_ylim(-90, 90)  # Adjust limits for second scale
    ax2.set_ylabel(r"$\alpha$ [deg]", color="tab:orange")
    ax2.tick_params(axis='y', labelcolor="tab:orange")

    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")

    plt.show()

"""
def plot_weight_distribution(w_vector1, w_vector2):
    barWidth = 0.5

    #      fus, wing,  tail
    atr = [w_vector1[0], w_vector1[1], w_vector1[2]]
    duuc = [w_vector2[0], w_vector2[1], w_vector2[2]]

    # Bar positions
    r = np.arange(len(atr))
    r2 = r + barWidth

    # Plotting
    fig, ax = plt.subplots(figsize=(6, 5))
    fig.canvas.manager.set_window_title("Weight per component")
    ax.bar(r, atr, color='tab:blue', width=barWidth, edgecolor='white', label='ATR72-600')
    ax.bar(r2, duuc, color='tab:orange', width=barWidth, edgecolor='white', label='DUUC')

    # Xticks
    ax.set_title("Weight comparison")
    ax.set_ylabel('component mass [kg]')
    ax.set_xticks(r + barWidth)
    ax.set_xticklabels(['Fuselage', 'Wing', 'Tail'])
    ax.legend()
    plt.show() """


def plot_weight_distribution(w_vector1, w_vector2):
    fus_duuc = w_vector1[0]
    fus_atr = w_vector2[0]

    wing_duuc = w_vector1[1]
    wing_atr = w_vector2[1]

    prop_nac_atr = w_vector2[2]
    htail_atr = w_vector2[3]
    vtail_atr = w_vector2[4]
    cv_atr = w_vector2[5]

    duct_duuc = w_vector1[2]
    pylon_duuc = w_vector1[3]
    support_duuc = w_vector1[4]
    cv_duuc = w_vector1[5]
    prop_duuc = w_vector1[6]  # also includes nacelle

    # Data
    labels = ['DUUC', 'ATR72-600']
    fuselage = np.array([fus_duuc, fus_atr])
    wing = np.array([wing_duuc, wing_atr])
    engine = np.array([0, prop_nac_atr])
    duct = np.array([duct_duuc, 0])
    pylon = np.array([pylon_duuc, 0])
    control_vanes = np.array([cv_duuc, cv_atr])
    support = np.array([support_duuc, 0])
    propeller = np.array([prop_duuc, 0])
    htail = np.array([0, htail_atr])
    vtail = np.array([0, vtail_atr])

    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.canvas.manager.set_window_title("Weight per component")

    fig.suptitle("Component Mass Comparison")

    # Plot 1: Fuselage
    axes[0].bar(labels, fuselage, label='Fuselage', color='tab:blue')
    axes[0].set_title('Fuselage')
    axes[0].legend(loc="lower right")

    # Plot 2: Wing + Engine (Stacked)
    axes[1].bar(labels, wing, label='Wing', color='tab:blue')
    axes[1].bar(labels, engine, label='Engine', bottom=wing, color='tab:orange')
    axes[1].set_title('Wing')
    axes[1].legend(loc="lower right")

    # Plot 3: Empennage (Stacked)
    axes[2].bar(labels, duct, label='Duct', color='tab:purple')
    axes[2].bar(labels, pylon, label='Pylon', bottom=duct, color='tab:olive')
    axes[2].bar(labels, control_vanes, label='Control Vanes', bottom=duct + pylon, color='tab:grey')
    axes[2].bar(labels, support, label='Support', bottom=duct + pylon + control_vanes, color='tab:red')
    axes[2].bar(labels, propeller, label='Propeller', bottom=duct + pylon + control_vanes + support, color='tab:green')
    axes[2].bar(labels, htail, label='Horizontal tail', bottom=duct + pylon + control_vanes + support + propeller, color='tab:orange')
    axes[2].bar(labels, vtail, label='Vertical tail', bottom=duct + pylon + control_vanes + support + propeller + htail, color='tab:blue')
    axes[2].set_title('Empennage')
    axes[2].legend(loc="lower right")
    axes[0].set_ylabel('Component mass [kg]')
    # Adjust layout and display
    plt.tight_layout()
    plt.show()

