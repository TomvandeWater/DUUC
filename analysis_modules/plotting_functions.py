import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import data.atr_reference as ref


def plot_inflow_properties2(v_input, a_input, station):
    v_input2 = [100, 150, 175, 150, 125, 150, 180, 100]
    a_input2 = [0, 10, 20, 30, 25, 10, 25, 0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.canvas.manager.set_window_title("Inflow parameters per component")

    # Load background image
    img = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\propulsive_empennage.png")

    # Plot for inflow velocity
    ax1.imshow(img, extent=[0.35, 6.00, 0.00, 200.00])
    ax1.plot(station, v_input, label=r"$V_{inflow}$ - power on", color="tab:blue", marker='o')
    ax1.plot(station, v_input2, label=r"$V_{inflow}$ - power off", color="tab:blue", linestyle="dotted", marker='o')
    ax1.set_xlim(0, 6)
    ax1.set_ylim(0, 200)
    ax1.xaxis.set_major_locator(ticker.NullLocator())
    ax1.set_aspect(aspect=0.02)
    ax1.set_ylabel(r'$V_{\infty}$ [m/s]')
    ax1.set_title(r'Inflow Velocity per Component')
    ax1.tick_params(axis='y')
    ax1.grid(True)
    ax1.legend(loc="upper left")

    # Plot for inflow angle
    ax2.imshow(img, extent=[0.35, 6.00, -100, 100])
    ax2.plot(station, a_input, label=r"$\alpha$ - power on", color="tab:orange", marker='o')
    ax2.plot(station, a_input2, label=r"$\alpha$ - power off", color="tab:orange", linestyle="dotted", marker='o')
    ax2.set_xlim(0, 6)
    ax2.set_ylim(-100, 100)
    ax2.xaxis.set_major_locator(ticker.NullLocator())
    ax2.set_aspect(aspect=0.02)
    ax2.set_ylabel(r"$\alpha$ [deg]")
    ax2.set_title(r'Inflow Angle per Component')
    ax2.tick_params(axis='y')
    ax2.grid(True)
    ax2.legend(loc="upper left")

    plt.tight_layout()
    plt.show()


def plot_weight_distribution(w_vector1, w_vector2):
    fus_duuc = w_vector1[0]
    fus_atr = w_vector2[0]

    wing_duuc = w_vector1[1]
    wing_atr = w_vector2[1]

    engine_atr = w_vector2[2]
    htail_atr = w_vector2[3]
    vtail_atr = w_vector2[4]
    cv_atr = w_vector2[5]
    nacelle_atr = w_vector2[6]
    fan_atr = w_vector2[8]

    duct_duuc = w_vector1[2]
    pylon_duuc = w_vector1[3]
    support_duuc = w_vector1[4]
    cv_duuc = w_vector1[5]
    engine_duuc = w_vector1[6]
    nacelle_duuc = w_vector1[7]
    fan_duuc = w_vector1[8]

    # Data
    labels = ['DUUC', 'ATR72-600']
    fuselage = np.array([fus_duuc, fus_atr])
    wing = np.array([wing_duuc, wing_atr])
    engine = np.array([0, engine_atr * 2])
    controls_wing = np.array([cv_duuc * 0.75, cv_atr * 0.75])
    duct = np.array([duct_duuc * 2, 0])
    pylon = np.array([pylon_duuc * 2, 0])
    control_vanes = np.array([cv_duuc * 0.25, cv_atr * 0.25])
    support = np.array([support_duuc * 2, 0])
    propeller = np.array([engine_duuc * 2, 0])
    htail = np.array([0, htail_atr])
    vtail = np.array([0, vtail_atr])
    nacelle_wing = np.array([0, nacelle_atr * 2])
    fan_wing = np.array([0, fan_atr * 2])
    nacelle_tail = np.array([nacelle_duuc * 2, 0])
    fan_tail = np.array([fan_duuc * 2, 0])

    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.canvas.manager.set_window_title("Weight per component")

    fig.suptitle("Component Mass Comparison")

    # Plot 1: Fuselage
    axes[0].bar(labels, fuselage, label='Fuselage', color='tab:blue')
    axes[0].set_title('Fuselage')
    axes[0].legend(loc="lower right")
    axes[0].set_ylim([0, 6000])

    # Plot 2: Wing + Engine (Stacked)
    axes[1].bar(labels, wing, label='Wing')
    axes[1].bar(labels, controls_wing, label='Controls', bottom=wing)
    axes[1].bar(labels, engine, label='Engine', bottom=wing + controls_wing)
    axes[1].bar(labels, nacelle_wing, label='Nacelle', bottom=wing + controls_wing + engine)
    axes[1].bar(labels, fan_wing, label='Fan', bottom=wing + controls_wing + engine + nacelle_wing)
    axes[1].set_ylim([0, 6000])

    axes[1].set_title('Wing')
    axes[1].legend(loc="lower right")

    # Plot 3: Empennage (Stacked)
    axes[2].bar(labels, control_vanes, label='Control Vanes')
    axes[2].bar(labels, duct, label='Duct', bottom=control_vanes)
    axes[2].bar(labels, pylon, label='Pylon', bottom=duct + control_vanes)
    axes[2].bar(labels, support, label='Support', bottom=duct + pylon + control_vanes)
    axes[2].bar(labels, propeller, label='Propeller', bottom=duct + pylon + control_vanes + support)
    axes[2].bar(labels, htail, label='Horizontal tail', bottom=duct + pylon + control_vanes + support + propeller)
    axes[2].bar(labels, vtail, label='Vertical tail', bottom=duct + pylon + control_vanes + support + propeller + htail)
    axes[2].bar(labels, nacelle_tail, label='Nacelle', bottom=duct + pylon + control_vanes + support + propeller
                                                                   + htail + vtail)
    axes[2].bar(labels, fan_tail, label='Fan', bottom=duct + pylon + control_vanes + support + propeller + htail
                                                           + vtail + nacelle_tail)
    axes[2].set_title('Empennage')
    axes[2].legend(loc="lower right")
    axes[2].set_ylim([0, 6000])
    axes[0].set_ylabel('Component mass [kg]')
    # Adjust layout and display
    plt.tight_layout()
    plt.show()


def weight_comp_table(w_vector, datatype: str):
    if datatype == "conventional":
        mtom = ref.MTOW
        ref_mat = [124.0, 178.0, 974.0, 766.0, 241.0, 0.0, 3045.0, 2323.0]

        components = [
            "Horizontal Tail", "Vertical Tail", "Controls", "Engine", "Nacelle",
            "Fan", "Wing", "Fuselage"
        ]

        computed_values = [
            w_vector[3], w_vector[4], w_vector[5], w_vector[2], w_vector[6],
            w_vector[8], w_vector[1], w_vector[0]
        ]

        data = []
        sum_perc = 0
        for i in range(len(components)):
            code_value = np.round(computed_values[i], 0)
            delta = np.round(code_value - ref_mat[i], 0)
            percent_delta = round(100 * delta / mtom, 4)
            data.append((components[i], ref_mat[i], code_value, delta, percent_delta))
            sum_perc += percent_delta

        header = "{:<20} | {:<20} | {:<20} | {:<20} | {:<10} |".format(
            "Component", "Reference Value [kg]", "Code Value [kg]", "Delta [kg]", "% Delta")
        separator = "-" * len(header)

        print("\n----- WEIGHT COMPARISON ATR72-600 WITH REFERENCE -----")
        print(separator)
        print(header)
        print(separator)

        for item, ref_val, code_val, delta, percent_delta in data:
            print("{:<20} | {:<20} | {:<20} | {:<20} | {:<10.4f} |".format(
                item, ref_val, code_val, delta, percent_delta))
        print(separator)
        print("{:<88} | {} %".format("Total deviation", np.round(sum_perc, 4)))
    if datatype == 'DUUC':
        ref_dungen = [265, 359, 0, 0, 81, 196, 655, 2569, 2685]
        ref_stavreva = [236.1, 115.4, 0, 0, 182, 219, 1315, 3310, 3507]
        mtom = ref.MTOW

        in_vector = [np.round(w_vector[2]), np.round(w_vector[3]), np.round(w_vector[5]), np.round(w_vector[4]),
                     np.round(w_vector[8]), np.round(w_vector[7]), np.round(w_vector[6]), np.round(w_vector[1]),
                     np.round(w_vector[0])]

        d_vector = [[np.round(in_vector[i] - ref_dungen[i], 0), np.round(in_vector[i] - ref_stavreva[i], 0)] for i in
                    range(len(ref_dungen))]
        mtom_vector = [[round(d_vector[i][0] / mtom * 100, 2), round(d_vector[i][1] / mtom * 100, 2)] for i in
                       range(len(ref_dungen))]

        components = ["Duct", "Pylon", "Controls", "Support", "Fan", "Nacelle", "Engine", "Wing", "Fuselage"]

        data = [
            (comp, ref_dungen[i], ref_stavreva[i], in_vector[i], d_vector[i][0], d_vector[i][1], mtom_vector[i][0],
             mtom_vector[i][1])
            for i, comp in enumerate(components)
        ]

        header = "{:<15} | {:<20} | {:<20} | {:<15} | {:<15} | {:<15} | {:<10} | {:<10} |".format(
            "Component", "Ref Dungen [kg]", "Ref Stavreva [kg]", "Code Value [kg]", "Delta Dungen", "Delta Stavreva",
            "% Dungen", "% Stavreva"
        )
        separator = "-" * len(header)

        print("\n----- WEIGHT COMPARISON DUUC WITH PREVIOUS WORK -----")
        print(separator)
        print(header)
        print(separator)

        sum_percd = 0
        sum_perc_str = 0

        for row in data:
            print("{:<15} | {:<20} | {:<20} | {:<15} | {:<15} | {:<15} | {:<10.2f} | {:<10.2f} |".format(*row))
            sum_percd += row[6]
            sum_perc_str += row[7]
        print(separator)
        print("{:<115} | {:<8} % | {} %".format("Total deviation", np.round(sum_percd, 4), np.round(sum_perc_str, 4)))
    else:
        return None


def interference_drag_range(cdi_results):
    cdi_results_matrix = np.array(cdi_results)

    plt.figure("Interference drag coefficients")
    plt.plot(cdi_results_matrix[:, 0], cdi_results_matrix[:, 7], label=r'Propeller')
    plt.plot(cdi_results_matrix[:, 0], cdi_results_matrix[:, 1], label=r'Control Vanes')
    plt.plot(cdi_results_matrix[:, 0], cdi_results_matrix[:, 4], label=r'Support')
    plt.plot(cdi_results_matrix[:, 0], cdi_results_matrix[:, 3], label=r'Pylon')
    plt.xlabel(r'$V_{\infty}$ [m/s]')
    plt.ylabel(r'$C_{d_{interference}}$ [-]')
    plt.title(r'Interference drag coefficients')
    plt.grid(True)
    plt.legend()
    plt.show()


def component_drag_wcdi_range(cd_results, cdi_results, cd0_results, j, component: str):
    # propeller j = 7
    # nacelle j = 6
    # cv  j =5
    # Duct = 2

    cdi_results_matrix = np.array(cdi_results)
    cd_results_matrix = np.array(cd_results)
    cd0_results_matrix = np.array(cd0_results)

    plt.figure(f"{component} drag coefficients")
    plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, j], label=r'Total drag coefficient')
    plt.plot(cd_results_matrix[:, 0], cdi_results_matrix[:, j], label=r'Induced drag coefficient')
    plt.plot(cd_results_matrix[:, 0], cd0_results_matrix[:, j], label=r'Zero lift drag coefficient')
    plt.xlabel(r'$v_{inf} [m/s]$')
    plt.ylabel(r'$C_{d} [-]$')
    plt.title(f'{component} drag coefficients')
    plt.legend()
    plt.show()


def component_drag_range(cd_results, cd0_results, j, component: str):
    # support j = 4
    # pylon j = 3

    cd_results_matrix = np.array(cd_results)
    cd0_results_matrix = np.array(cd0_results)

    plt.figure(f"{component} drag coefficients")
    plt.plot(cd_results_matrix[:, 0], cd_results_matrix[:, j], label=r'Total drag coefficient')

    plt.plot(cd_results_matrix[:, 0], cd0_results_matrix[:, j], label=r'Zero lift drag coefficient')
    plt.xlabel(r'$v_{inf} [m/s]$')
    plt.ylabel(r'$C_{d} [-]$')
    plt.title(f'{component} drag coefficients')
    plt.legend()


def propulsive_efficiency_plot(eta_prop_results):
    eta_results_matrix = np.array(eta_prop_results)

    plt.figure('Propulsive efficiency')
    plt.plot(eta_results_matrix[:, 0], eta_results_matrix[:, 1], label=r'DUUC')
    plt.plot(eta_results_matrix[:, 0], eta_results_matrix[:, 2], label=r'ATR72-600')
    plt.xlabel(r'J [-]')
    plt.ylabel(r'$\eta_{prop}$ [%]')
    plt.title(r'Propulsive efficiency')
    plt.legend()
    plt.grid(True)
