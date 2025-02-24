import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker


def plot_inflow_velocity(v_input, a_input, station):
    #station = [0, 190, 420, 650, 1100, 1300]

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


def plot_inflow_angle(a_input, station):
    plt.figure("Inflow angle per component")
    plt.plot(station, a_input)
    plt.xlabel(r'Station numbers')
    plt.ylabel(r'$\alpha$ [deg]')
    plt.title(r'Inflow angle per component')
    # plt.legend()
    plt.grid(True)
    plt.show()


def plot_weight_distribution(w_vector1, w_vector2):
    barWidth = 0.5

    #      fus, wing, prop, tail
    atr = [w_vector1[0], w_vector1[1], w_vector1[2], w_vector1[3]]
    duuc = [w_vector2[0], w_vector2[1], w_vector2[2], w_vector2[3]]

    # Bar positions
    r = np.arange(len(atr))
    r2 = r + barWidth

    # Plotting
    fig, ax = plt.subplots(dpi=300)
    ax.bar(r, atr, color='#7f6d5f', width=barWidth, edgecolor='white', label='ATR72-600')
    ax.bar(r2, duuc, color='#557f2d', width=barWidth, edgecolor='white', label='DUUC')

    # Xticks
    ax.set_xlabel('group', fontweight='bold')
    ax.set_xticks(r + barWidth)
    ax.set_xticklabels(['Fuselage', 'Wing', 'Propulsion', 'Tail'])
    ax.legend()
    plt.show()
