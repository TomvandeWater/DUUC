import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np


# Function to generate plots for each component
def create_plot(ax, component_name, plot_type, alpha, CL, CD, CM, reference_alpha=None, reference_CL=None,
                reference_CD=None,
                reference_CM=None, reference_alpha2=None, reference_CL2=None, reference_CD2=None, reference_CM2=None):
    # Plot the main results
    if plot_type == 'CL_vs_Alpha':
        ax.plot(alpha, CL, label="Prediction Model", color='blue')
        if reference_alpha is not None and reference_CL is not None and component_name != "Duct":
            ax.plot(reference_alpha, reference_CL, label="XFoil", color='blue', linestyle='dashed')

        if reference_alpha is not None and reference_CL is not None and component_name == "Duct":
            ax.plot(reference_alpha, reference_CL, label="Literature reference 2", color='green', linestyle='dashed')

        if reference_alpha2 is not None and reference_CL2 is not None:
            ax.plot(reference_alpha2, reference_CL2, label="Literature reference", color='blue', marker='o')
        ax.set_title('CL vs. Alpha')
        ax.set_xlabel(r'$\alpha$')
        ax.set_ylabel(r'$C_{L}$')
        ax.legend(loc='upper left')

    elif plot_type == 'CD_vs_Alpha':
        ax.plot(alpha, CD, label="Prediction Model", color='red')
        if reference_alpha is not None and reference_CD is not None and component_name != "Duct":
            ax.plot(reference_alpha, reference_CD, label="XFoil", color='red', linestyle='dashed')

        if reference_alpha is not None and reference_CD is not None and component_name == "Duct":
            ax.plot(reference_alpha, reference_CL, label="Literature reference 2", color='green', linestyle='dashed')

        if reference_alpha2 is not None and reference_CD2 is not None:
            ax.plot(reference_alpha2, reference_CD2, label="Literature reference", color='red', marker='o')
        ax.set_title('CD vs. Alpha')
        ax.set_xlabel(r'$\alpha$')
        ax.set_ylabel(r'$C_{D}$')
        ax.legend(loc='upper left')

    elif plot_type == 'CL_vs_CD':
        ax.plot(CD, CL, label="Prediction Model", color='green')
        if reference_CD is not None and reference_CL is not None and component_name != "Duct":
            ax.plot(reference_CD, reference_CL, label="XFoil", color='green', linestyle='dashed')

        if reference_CD is not None and reference_CL is not None and component_name == "Duct":
            ax.plot(reference_CD, reference_CL, label="Literature reference 2", color='green', linestyle='dashed')

        if reference_CD2 is not None and reference_CL2 is not None:
            ax.plot(reference_CD2, reference_CL2, label="Literature reference", color='green', marker='o')
        ax.set_title('CL vs. CD')
        ax.set_xlabel(r'$C_{D}$')
        ax.set_ylabel(r'$C_{L}$')
        ax.legend(loc='upper left')

    elif plot_type == 'CM_vs_Alpha':
        ax.plot(alpha, CM, label="Prediction Model", color='purple')
        if reference_alpha is not None and reference_CM is not None and component_name != "Duct":
            ax.plot(reference_alpha, reference_CM, label="XFoil", color='purple', linestyle='dashed')

        if reference_alpha is not None and reference_CM is not None and component_name == "Duct":
            ax.plot(reference_alpha, reference_CM, label="Literature reference 2", color='green', linestyle='dashed')

        if reference_alpha2 is not None and reference_CM2 is not None:
            ax.plot(reference_alpha2, reference_CM2, label="Literature reference", color='purple', marker='o')
        ax.set_title(r'CM vs. $\alpha$')
        ax.set_xlabel(r'$\alpha$')
        ax.set_ylabel(r'$C_{M}$')
        ax.legend(loc='upper left')

    elif plot_type == 'CD_vs_CL_squared':
        ax.plot([cl ** 2 for cl in CL], CD, label=r'Model', color="tab:blue")
        ax.set_title(r'$C_D$ vs. $C_L^2$')
        ax.set_xlabel(r'$C_L^2$')
        ax.set_ylabel(r'$C_D$')
        ax.legend(loc='upper left')


# Function to add a set of 4 plots and a table into a tab for a component
def add_plots_to_tab(tab, component_name, results_component, reference_results=None, reference_results2=None):
    # Unpack results_component for plotting
    alpha = []
    CL = []
    CD = []
    CM = []

    for result in results_component:
        a, cl, cd, cm = result  # Unpack the result into individual values
        alpha.append(a)
        CL.append(cl)
        CD.append(cd)
        CM.append(cm)

    # If reference results are provided, unpack them as well
    reference_alpha = []
    reference_CL = []
    reference_CD = []
    reference_CM = []

    if reference_results:
        for ref_result in reference_results:
            ref_a, ref_cl, ref_cd, ref_cm = ref_result
            reference_alpha.append(ref_a)
            reference_CL.append(ref_cl)
            reference_CD.append(ref_cd)
            reference_CM.append(ref_cm)

    reference_alpha2 = []
    reference_CL2 = []
    reference_CD2 = []
    reference_CM2 = []

    if reference_results2:
        for ref_result2 in reference_results2:
            ref_a, ref_cl, ref_cd, ref_cm = ref_result2
            reference_alpha2.append(ref_a)
            reference_CL2.append(ref_cl)
            reference_CD2.append(ref_cd)
            reference_CM2.append(ref_cm)

    # Convert lists to numpy arrays for plotting
    alpha = np.array(alpha)
    CL = np.array(CL)
    CD = np.array(CD)
    CM = np.array(CM)
    reference_alpha = np.array(reference_alpha) if reference_alpha else None
    reference_CL = np.array(reference_CL) if reference_CL else None
    reference_CD = np.array(reference_CD) if reference_CD else None
    reference_CM = np.array(reference_CM) if reference_CM else None

    reference_alpha2 = np.array(reference_alpha2) if reference_alpha2 else None
    reference_CL2 = np.array(reference_CL2) if reference_CL2 else None
    reference_CD2 = np.array(reference_CD2) if reference_CD2 else None
    reference_CM2 = np.array(reference_CM2) if reference_CM2 else None

    # Create the figure and axes for 2x3 layout
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))  # 2x3 grid for 6 plots

    # Flatten the axs array for easier iteration
    axs = axs.flatten()

    # Create and add each plot for the component
    plot_types = ['CL_vs_Alpha', 'CD_vs_Alpha', 'CL_vs_CD', 'CM_vs_Alpha']
    for i, plot_type in enumerate(plot_types):
        create_plot(axs[i], component_name, plot_type, alpha, CL, CD, CM, reference_alpha, reference_CL, reference_CD,
                    reference_CM,
                    reference_alpha2, reference_CL2, reference_CD2, reference_CM2)

    # Create the new plot for CD vs CL^2 and add it to axs[4] (center plot)
    create_plot(axs[4], component_name, 'CD_vs_CL_squared', alpha, CL, CD, CM, reference_alpha, reference_CL,
                reference_CD, reference_CM,
                reference_alpha2, reference_CL2, reference_CD2, reference_CM2)

    # Add the table in the sixth plot (axs[5])
    table_data = [
        ['Parameter', 'Value'],
        ['Alpha', f'{alpha.mean():.2f}'],
        ['CL', f'{CL.mean():.2f}'],
        ['CD', f'{CD.mean():.2f}'],
        ['CM', f'{CM.mean():.2f}']
    ]

    # Place the table inside the subplot at axs[5] (bottom-right of the 2x3 grid)
    table = axs[5].table(cellText=table_data, loc='center', cellLoc='center', colWidths=[0.3, 0.3])
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.2)
    axs[5].axis('off')  # Hide the axis

    # Adjust layout to avoid title overlap
    plt.tight_layout()  # Automatically adjust spacing

    # Embed the matplotlib figure into the Tkinter tab
    canvas = FigureCanvasTkAgg(fig, master=tab)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


# Function to create the Tkinter GUI
def create_gui(results_duct, results_pylon, results_nacelle, results_support, results_control_vanes,
               reference_results_duct=None, reference_results_pylon=None, reference_results_nacelle=None,
               reference_results_support=None, reference_results_control_vanes=None,
               reference_results_duct2=None, reference_results_pylon2=None, reference_results_nacelle2=None,
               reference_results_support2=None, reference_results_control_vanes2=None
               ):
    root = tk.Tk()
    root.title("Aerodynamic Component Plots")

    # Create a notebook (tab container)
    notebook = ttk.Notebook(root)
    notebook.pack(fill=tk.BOTH, expand=True)

    # Component names for the tabs
    component_names = ['Duct', 'Pylon', 'Nacelle', 'Support', 'Control Vanes']

    # Add tabs to the notebook with custom names
    results_dict = {
        'Duct': results_duct,
        'Pylon': results_pylon,
        'Nacelle': results_nacelle,
        'Support': results_support,
        'Control Vanes': results_control_vanes
    }

    reference_results_dict = {
        'Duct': reference_results_duct,
        'Pylon': reference_results_pylon,
        'Nacelle': reference_results_nacelle,
        'Support': reference_results_support,
        'Control Vanes': reference_results_control_vanes
    }
    reference_results_dict2 = {
        'Duct': reference_results_duct2,
        'Pylon': reference_results_pylon2,
        'Nacelle': reference_results_nacelle2,
        'Support': reference_results_support2,
        'Control Vanes': reference_results_control_vanes2
    }

    for component_name in component_names:
        tab = ttk.Frame(notebook)
        notebook.add(tab, text=component_name)  # Set the tab name

        # Retrieve the results for the component from the passed data
        results_component = results_dict.get(component_name, [])
        reference_results = reference_results_dict.get(component_name, None)
        reference_results2 = reference_results_dict2.get(component_name, None)

        # Add plots for each component to its respective tab
        add_plots_to_tab(tab, component_name, results_component, reference_results, reference_results2)

    root.mainloop()



