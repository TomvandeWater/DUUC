import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import data.atr_reference as ref
from analysis_modules.aerodynamic import cl_cd_bucket_polar


def plot_inflow_properties2(v_input, a_input, station):
    v_input2 = [127.0, 127.0, 78, 78, 78, 63, 96, 127.0]
    a_input2 = [0, 0, 0, 0, 0, 0, 0, 0]

    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    # fig.canvas.manager.set_window_title("Inflow parameters per component")

    # Load background image
    img = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\propulsive_empennage.png")

    # Plot for inflow velocity
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    fig1.canvas.manager.set_window_title("Inflow Velocity per Component")
    ax1.imshow(img, extent=[0.35, 6.00, 0.00, 200.00], aspect='auto')
    ax1.plot(station, v_input, label=r"$V_{inflow}$ - power on", color="tab:blue", marker='o')
    ax1.plot(station, v_input2, label=r"$V_{inflow}$ - power off", color="tab:blue", linestyle="dotted", marker='o')
    ax1.set_xlim(0, 6)
    ax1.set_ylim(0, 200)
    ax1.xaxis.set_major_locator(ticker.NullLocator())
    #ax1.set_aspect(aspect=0.02)
    ax1.set_ylabel(r'$V_{\infty}$ [m/s]')
    ax1.set_title(r'Inflow Velocity per Component')
    ax1.tick_params(axis='y')
    ax1.grid(True)
    ax1.legend(loc="upper left")
    #plt.tight_layout

    # Plot for inflow angle
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    fig2.canvas.manager.set_window_title("Inflow Angle per Component")
    ax2.imshow(img, extent=[0.35, 6.00, -100, 100], aspect='auto')
    ax2.plot(station, a_input, label=r"$\alpha$ - power on", color="tab:orange", marker='o')
    ax2.plot(station, a_input2, label=r"$\alpha$ - power off", color="tab:orange", linestyle="dotted", marker='o')
    ax2.set_xlim(0, 6)
    ax2.set_ylim(-100, 100)
    ax2.xaxis.set_major_locator(ticker.NullLocator())
    #ax2.set_aspect(aspect=0.02)
    ax2.set_ylabel(r"$\alpha$ [deg]")
    ax2.set_title(r'Inflow Angle per Component')
    ax2.tick_params(axis='y')
    ax2.grid(True)
    ax2.legend(loc="upper left")

    #plt.tight_layout()
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


def cd0_drag_comparison(cd0_vector1, cd0_vector2):
    cd0_vector1 = np.array(cd0_vector1)
    cd0_ref1 = np.array([1.6e-3, 8.347e-4, 2.615e-3, 0.014, 8.053e-3])
    cd0_vector2 = np.array(cd0_vector2)
    cd0_ref2 = np.array([0, 0, 0, 0, 0, 0, 0])
    cd_totals = np.array([sum(cd0_vector1), sum(cd0_vector2)])
    cd_totals_ref = np.array([sum(cd0_ref1), sum(cd0_ref2)])

    # Labels for the x-axis
    label1 = ['Nacelle', 'Horizontal Tail', 'Vertical Tail', 'Wing', 'Fuselage']
    label2 = ['Duct', 'Pylon', 'Nacelle', 'Support', 'Control', 'Wing', 'Fuselage']
    label3 = ['ATR72-600', 'DUUC']

    # Create the figure and subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    fig.canvas.manager.set_window_title("Zero Lift Drag comparison")
    fig.suptitle("Zero Lift Drag per component")

    x = np.arange(len(label1))  # the label locations
    width = 0.35  # the width of the bars

    ax1.bar(x - width / 2, cd0_vector1, width, label='Model', color="tab:orange")
    ax1.bar(x + width / 2, cd0_ref1, width, label='Reference', color="tab:green")  # added ref bar and label.

    ax1.set_title('ATR72-600')
    ax1.set_ylabel('$C_{D0}$ [-]')
    ax1.set_ylim([0, 0.06])
    ax1.set_xticks(x)
    ax1.set_xticklabels(label1)
    ax1.legend()

    # Plot the second subplot
    x = np.arange(len(label2))
    ax2.bar(x - width / 2, cd0_vector2, width, label='Model', color="tab:blue")
    ax2.bar(x + width / 2, cd0_ref2, width, label='Reference', color="tab:green")  # added ref bar and label.
    ax2.set_title('DUUC')
    ax2.set_ylabel('$C_{D0}$ [-]')
    ax2.set_ylim([0, 0.06])
    ax2.set_xticks(x)
    ax2.set_xticklabels(label2)
    ax2.legend()

    colors = ['tab:orange', 'tab:blue']
    x = np.arange(len(label3))
    ax3.bar(x - width / 2, cd_totals, width, label='Model', color=colors)
    ax3.bar(x + width / 2, cd_totals_ref, width, label='Reference', color="tab:green")  # added ref bar and label.
    ax3.set_title('Aircraft sum')
    ax3.set_ylabel('$C_{D0}$ [-]')
    ax3.set_xticks(x)
    ax3.set_xticklabels(label3)
    ax3.legend()

    # Adjust layout to prevent overlapping titles/labels
    plt.tight_layout()

    # Show the plot
    plt.show()


def cd0_drag_empennage(cd0_vector1, cd0_vector2):
    cd0_vector1 = np.array(cd0_vector1)
    cd0_ref1 = np.array([8.347e-4, 2.065e-3])
    cd0_vector2 = np.array(cd0_vector2)
    cd0_ref2 = np.array([49.2e-4, 2.8e-4, 0, 0])
    cd_totals = np.array([sum(cd0_vector1), sum(cd0_vector2)])
    cd_totals_ref = np.array([sum(cd0_ref1), sum(cd0_ref2)])

    # Labels for the x-axis
    label1 = ['Hor. Tail', 'Vert. Tail']
    label2 = ['Duct', 'Pylon', 'Support', 'Control']
    label3 = ['ATR72-600', 'DUUC']

    # Combine labels for both ATR72-600 and DUUC
    combined_labels = label1 + label2
    combined_cd0_vector = np.concatenate([cd0_vector1, cd0_vector2])
    combined_cd0_ref = np.concatenate([cd0_ref1, cd0_ref2])

    # Create the figure and subplots
    # fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(10, 6))
    # fig.canvas.manager.set_window_title("Zero Lift Drag comparison")
    # fig.suptitle("Zero Lift Drag per component")

    width = 0.35  # the width of the bars

    fig1, ax1 = plt.subplots()
    fig1.canvas.manager.set_window_title("Zero lift drag per component")
    # Plot combined components of ATR72-600 and DUUC in the first subplot
    x = np.arange(len(combined_labels))  # the label locations
    ax1.bar(x - width / 2, combined_cd0_vector, width, label='Model', color="tab:blue")
    ax1.bar(x + width / 2, combined_cd0_ref, width, label='Reference', color="tab:green")

    ax1.set_title('Comparison of Empennage components')
    ax1.set_ylabel('$C_{D0}$ [-]')
    ax1.set_ylim([0, 0.01])
    ax1.set_xticks(x)
    ax1.set_xticklabels(combined_labels)
    ax1.legend(loc='upper left')

    fig2, ax3 = plt.subplots()
    fig2.canvas.manager.set_window_title("Zero lift drag empennage sum")
    colors = ['tab:blue', 'tab:blue']
    x = np.arange(len(label3))
    ax3.bar(x - width / 2, cd_totals, width, label='Model', color=colors)
    ax3.bar(x + width / 2, cd_totals_ref, width, label='Reference', color="tab:green")  # added ref bar and label.
    ax3.set_title('Empennage totals')
    ax3.set_ylabel('$C_{D0}$ [-]')
    ax3.set_xticks(x)
    ax3.set_xticklabels(label3)
    ax3.legend(loc='upper left')

    # Adjust layout to prevent overlapping titles/labels
    #plt.tight_layout()

    # Show the plot
    plt.show()



def cl_cd_bucket(cd01, cdi1, alpha):
    cl = np.linspace(-1.5, 3.5, 100)
    cd1 = []
    cd2 = []
    cd3 = []
    for i in range(len(cl)):
        cd1.append(cl_cd_bucket_polar(cd01, cdi1, cl[i]))
        cd2.append(cl_cd_bucket_polar(0.027403, 0.034, cl[i]))
        cd3.append(cl_cd_bucket_polar(0.050, 0.031, cl[i]))

    plt.figure(f'CL-CD bucket ({alpha})')
    plt.plot(cd1, cl, label=r'Model: ATR72-600', color="tab:orange")
    plt.plot(cd2, cl, label=r'Reference cruise', color="tab:green")
    plt.plot(cd3, cl, label='Reference take-off', color="tab:green", linestyle="--")
    plt.xlabel(r'$C_{D}$ [-]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.title(f'CL - CD Bucket ({alpha})')
    plt.legend()
    plt.grid(True)
    plt.show()

def cd_interference_drag_comparison(cd_vector1, cd_vector2, alpha, v_inf):
    cd0_vector1 = np.array(cd_vector1)
    cd0_ref1 = np.array([0, 0, 0])
    cd0_vector2 = np.array(cd_vector2)
    cd0_ref2 = np.array([0, 0, 0])
    cd_totals = np.array([sum(cd0_vector1), sum(cd0_vector2)])
    cd_totals_ref = np.array([sum(cd0_ref1), sum(cd0_ref2)])

    # Labels for the x-axis
    label1 = ['Nacelle', 'Horizontal Tail', 'Vertical Tail']
    label2 = ['Pylon', 'Support', 'Control']
    label3 = ['ATR72-600', 'DUUC']

    # Create the figure and subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    fig.canvas.manager.set_window_title("Interference Drag comparison")
    fig.suptitle(f"Interference per component (a = {alpha}, v = {v_inf} )")

    x = np.arange(len(label1))  # the label locations
    width = 0.35  # the width of the bars

    ax1.bar(x - width / 2, cd0_vector1, width, label='Model', color="tab:orange")
    ax1.bar(x + width / 2, cd0_ref1, width, label='Reference', color="tab:green")  # added ref bar and label.

    ax1.set_title('ATR72-600')
    ax1.set_ylabel('$C_{D-interference}$ [-]')
    ax1.set_ylim([0, 0.01])
    ax1.set_xticks(x)
    ax1.set_xticklabels(label1)
    ax1.legend()

    # Plot the second subplot
    x = np.arange(len(label2))
    ax2.bar(x - width / 2, cd0_vector2, width, label='Model', color="tab:blue")
    ax2.bar(x + width / 2, cd0_ref2, width, label='Reference', color="tab:green")  # added ref bar and label.
    ax2.set_title('DUUC')
    ax2.set_ylabel('$C_{D-interference}$ [-]')
    ax2.set_ylim([0, 0.01])
    ax2.set_xticks(x)
    ax2.set_xticklabels(label2)
    ax2.legend()

    colors = ['tab:orange', 'tab:blue']
    x = np.arange(len(label3))
    ax3.bar(x - width / 2, cd_totals, width, label='Model', color=colors)
    ax3.bar(x + width / 2, cd_totals_ref, width, label='Reference', color="tab:green")  # added ref bar and label.
    ax3.set_title('Aircraft sum')
    ax3.set_ylabel('$C_{D-interference}$ [-]')
    ax3.set_xticks(x)
    ax3.set_xticklabels(label3)
    ax3.legend()

    # Adjust layout to prevent overlapping titles/labels
    plt.tight_layout()

    # Show the plot
    plt.show()


def cl_comparison(cl_vector1, cl_vector2, alpha):
    cl_vector1 = np.array(cl_vector1)
    cl_ref1 = np.array([1.6e-3, 8.347e-4, 1.315e-3])
    cl_vector2 = np.array(cl_vector2)
    cl_ref2 = np.array([0, 0, 0, 0, 0, 0])
    cl_totals = np.array([sum(cl_vector1), sum(cl_vector2)])
    cl_totals_ref = np.array([sum(cl_ref1), sum(cl_ref2)])

    # Labels for the x-axis
    label1 = ['Horizontal Tail', 'Wing', 'Fuselage']
    label2 = ['Duct', 'Pylon', 'Support', 'Control', 'Wing', 'Fuselage']
    label3 = ['ATR72-600', 'DUUC']

    # Create the figure and subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    fig.canvas.manager.set_window_title(r"Cl comparison ")
    fig.suptitle(r"$C_L$ per component ($\alpha$ = {})".format(alpha))

    x = np.arange(len(label1))  # the label locations
    width = 0.35  # the width of the bars

    ax1.bar(x - width / 2, cl_vector1, width, label='Model', color="tab:orange")
    ax1.bar(x + width / 2, cl_ref1, width, label='Reference', color="tab:green")  # added ref bar and label.

    ax1.set_title('ATR72-600')
    ax1.set_ylabel('$C_{L}$ [-]')
    # ax1.set_ylim([0, 0.06])
    ax1.set_xticks(x)
    ax1.set_xticklabels(label1)
    ax1.legend()

    # Plot the second subplot
    x = np.arange(len(label2))
    ax2.bar(x - width / 2, cl_vector2, width, label='Model', color="tab:blue")
    ax2.bar(x + width / 2, cl_ref2, width, label='Reference', color="tab:green")  # added ref bar and label.
    ax2.set_title('DUUC')
    ax2.set_ylabel('$C_{L}$ [-]')
    # ax2.set_ylim([0, 0.06])
    ax2.set_xticks(x)
    ax2.set_xticklabels(label2)
    ax2.legend()

    colors = ['tab:orange', 'tab:blue']
    x = np.arange(len(label3))
    ax3.bar(x - width / 2, cl_totals, width, label='Model', color=colors)
    ax3.bar(x + width / 2, cl_totals_ref, width, label='Reference', color="tab:green")  # added ref bar and label.
    ax3.set_title('Aircraft sum')
    ax3.set_ylabel('$C_{L}$ [-]')
    ax3.set_xticks(x)
    ax3.set_xticklabels(label3)
    ax3.legend()

    # Adjust layout to prevent overlapping titles/labels
    plt.tight_layout()

    # Show the plot
    plt.show()


def print_cg_mass(part, mass, cg, mcg, m_sum, mcg_cum, cg_tot, group, aircraft):
    separator = "------------------------------------------------------------------------------------------------"
    header = "| Mass {} group  | Mass [kg]              | CG [m]                | Mass * CG [kg*m]     |"

    print("\n")
    print(f"---------------------------------------{aircraft}-----------------------------------------------")
    print(separator)
    print(header.format(group))
    print(separator)
    for i in range(len(mass)):
        print("| {:<20} | {:<22} | {:<21} | {:<20} |".format(part[i], mass[i], np.round(cg[i], 2),
                                                             np.round(mcg[i], 1)))
    print(separator)
    print("| sum mass             | {:<22} |                       | {:<20} |".format(m_sum, mcg_cum))
    print(separator)
    print("|                      |                        |     estimation cg     | {:<16} [m] |".format(np.round(
        cg_tot, 2)))
    print(separator)


def xplot(a1, b1, a2, xcg, static_margin, cmac, x_ac, aircraft):

    x1 = np.linspace(-1, 0, 101)
    x2 = np.linspace(-1, 2, 101)

    stm = static_margin / 100
    x_cg_ac = (x1 + x_ac)

    y1 = (a1 * x_cg_ac * cmac) + b1

    # Compute stability requirement
    y2 = a2 * x2
    y3 = y2 + stm  # Apply static margin shift

    # Find intersection of blue dotted line (CG location) with green dashed line (static margin line)
    y_intersection = a2 * xcg + stm
    x_intersect = (y_intersection - b1) / (a1 * cmac) - x_ac

    print(r"S_{}/S_{} = {}".format("H", "W", np.round(y_intersection, 4)))
    print(f"new S_h = {y_intersection * ref.s_w}")
    print(f"cg range: {xcg - x_intersect}")

    fig, axs = plt.subplots(2, 2, figsize=(14, 7))
    ax1, ax2, ax3 = axs[0, 0], axs[0, 1], axs[1, 0]

    fig.canvas.manager.set_window_title(f'X-plot - {aircraft}')
    fig.suptitle(f'X-plot - {aircraft}')

    ax1.plot(x1, y1, label=r'Control', color="tab:orange")
    ax1.plot(x2, y3, label=r'Static margin', color="tab:green", linestyle="--")
    ax1.plot(x2, y2, label=r'Stability', color="tab:green")
    ax1.vlines(xcg, 0, y_intersection, color="tab:blue", linestyle="dotted")
    ax1.scatter(xcg, 0, color='tab:blue', zorder=5, label="CG")
    ax1.hlines(y_intersection, x_intersect, xcg, color='tab:red', linestyle='dotted', label='CG range')
    ax1.scatter(xcg, y_intersection, color='tab:red', zorder=3)  # Start marker
    ax1.scatter(x_intersect, y_intersection, color='tab:red', zorder=3)
    ax1.set_xlabel(r'Center of Gravity [m]')
    ax1.set_ylabel(r'$S_{H}/S_{W}$ [-]')
    ax1.set_title("Model")
    ax1.set_ylim([-0.6, 0.8])
    ax1.legend()
    ax1.grid(True)

    a1_ref = -0.4487
    b1_ref = 0.20768
    a2_ref = 0.305
    x_ac_ref = 0.25 * 2.2345
    xcg_ref = 0.3086
    cmac_ref = 2.2345

    x1_ref = np.linspace(-1, 0, 101)
    x2_ref = np.linspace(-1, 2, 101)

    x_cg_ac_ref = (x1_ref + x_ac_ref)

    y1_ref = (a1_ref * x_cg_ac_ref * cmac_ref) + b1_ref

    # Compute stability requirement
    y2_ref = a2_ref * x2_ref
    y3_ref = y2_ref + 0.05  # Apply static margin shift

    # Find intersection of blue dotted line (CG location) with green dashed line (static margin line)
    y_intersection_ref = a2_ref * xcg_ref + 0.05
    x_intersect_ref = (y_intersection_ref - b1_ref) / (a1_ref * cmac_ref) - x_ac_ref

    print(r"S_{}/S_{} reference = {}".format("H", "W", np.round(y_intersection_ref, 4)))
    print(f"new S_h = {y_intersection_ref * ref.s_w}")
    print(f"cg range reference: {xcg_ref - x_intersect_ref}")

    ax2.plot(x1_ref, y1_ref, label=r'Control', color="tab:orange")
    ax2.plot(x2_ref, y3_ref, label=r'Static margin', color="tab:green", linestyle="--")
    ax2.plot(x2_ref, y2_ref, label=r'Stability', color="tab:green")
    ax2.vlines(xcg_ref, 0, y_intersection_ref, color="tab:blue", linestyle="dotted")
    ax2.scatter(xcg_ref, 0, color='tab:blue', zorder=5, label="CG")
    ax2.hlines(y_intersection_ref, x_intersect_ref, xcg_ref, color='tab:red', linestyle='dotted', label='CG range')
    ax2.scatter(xcg_ref, y_intersection_ref, color='tab:red', zorder=3)  # Start marker
    ax2.scatter(x_intersect_ref, y_intersection_ref, color='tab:red', zorder=3)
    ax2.set_xlabel(r'Center of Gravity [m]')
    ax2.set_ylabel(r'$S_{H}/S_{W}$ [-]')
    ax2.set_ylim([-0.6, 0.8])
    ax2.set_title("Reference data")
    ax2.legend()
    ax2.grid(True)

    # Third subplot: Table of key values
    table_data = [
        ['a1', np.round(a1, 4), np.round(a1_ref, 4), np.round(((a1 - a1_ref) / a1_ref) * 100, 2)],
        ['b1', np.round(b1, 4), np.round(b1_ref, 4), np.round(((b1 - b1_ref) / b1_ref) * 100, 2)],
        ['a2', np.round(a2, 4), np.round(a2_ref, 4), np.round(((a2 - a2_ref) / a2_ref) * 100, 2)],
        [r'$x_{cg}$', np.round(xcg, 4), np.round(xcg_ref, 4), np.round(((xcg - xcg_ref) / xcg_ref) * 100, 2)],
        [r'$S_{H}/S_{W}$', np.round(y_intersection, 4), np.round(y_intersection_ref, 4), np.round(((y_intersection - y_intersection_ref) / y_intersection_ref) * 100, 2)],
    ]

    col_labels = ['Parameter', 'Value model', 'Value reference', '% difference']
    ax3.axis('tight')
    ax3.axis('off')
    ax3.table(cellText=table_data, colLabels=col_labels, loc='center', cellLoc='center')
    axs[1, 1].axis('off')
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plt.show()

# xplot(-0.4887, 0.20768, 0.305, 0.3, 5, 2.2345, (0.25 * 2.2345), 'test')
