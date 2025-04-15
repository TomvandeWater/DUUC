import matplotlib.ticker as ticker
import os
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt6.QtWidgets import QWidget, QGridLayout, QTabWidget, QVBoxLayout
from analysis_modules.factors import *
from PyQt6.QtGui import QIcon
from data.read_data import get_polar
from analysis_modules.factors import area_ratio
from table_config import style_table
import data.experiment_reference_5annular_airfoil as ref5r
import matplotlib.gridspec as gridspec


class PlotManager(QWidget):
    def __init__(self, gui):
        super().__init__()
        self.gui = gui
        self.parameters = gui.parameters
        self.calculation_results = gui.calculation_results
        # print(f"calculation results: {self.calculation_results}")
        self.plots = {}  # Dictionary to store plots for each component

        # Define all plot components
        self.components = [
            "Inflow", "Geometry", "Duct", "Pylon", "Nacelle", "Support", "Control", "Empennage",
            "CD0", "Weight", "X_cog", "Vtail", "Htail", "requirements",
        ]

        # Initialize the layout
        self.layout = QGridLayout(self)

        # Create a tab widget
        self.tab_widget = QTabWidget()
        self.layout.addWidget(self.tab_widget, 0, 0)
        self.output_filepath = r'C:\Users\tomva\pythonProject\DUUC\data\plot_outputs'

        # Initialize tabs
        self.init_tabs()
        self.placeholder_icon = QIcon(r"C:\Users\tomva\pythonProject\DUUC\data\images\flame.ico")
        self.tab_widget.currentChanged.connect(self.on_tab_clicked)  # Connect to tab click

    def add_placeholder_icons(self, component):
        """Adds a placeholder icon to the center of each plot in the given component's tab."""
        for fig, ax, canvas in self.plots[component]:
            ax.clear()
            img = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\flame.ico")  # Load the placeholder image
            ax.imshow(img, aspect='auto', extent=[-1, 1, -1, 1])  # Set fixed aspect and limits
            ax.set_xlim(-1, 1)  # Fix x-axis limits
            ax.set_ylim(-1, 1)  # Fix y-axis limits
            ax.axis('off')  # Hide axes
            canvas.draw()

    def on_tab_clicked(self, index):
        """Handles the tab click event."""
        component = self.components[index]  # Get the component name from the tab index
        self.update_plots(component)  # Update plots on tab click

    def init_tabs(self):
        """Initialize tabs for all components."""
        self.cleanup()  # Clean up old plots before initializing new ones

        for component in self.components:
            tab = QWidget()
            tab_layout = QGridLayout(tab)

            if component == "X_cog":
                nested_tab_widget = QTabWidget(tab)
                tab_layout.addWidget(nested_tab_widget, 0, 0)
                x_cog_tabs = []
                tab_names = ["Center of Gravity", "Wing group", "Fuselage group", "Xcg with Xpe"]
                for j in range(4):
                    sub_tab = QWidget()
                    sub_tab_layout = QVBoxLayout(sub_tab)

                    fig, ax = plt.subplots(figsize=(12, 6) if j == 0 else (4, 3))  # Wider figure for first plot
                    canvas = FigureCanvas(fig)
                    sub_tab_layout.addWidget(canvas)
                    x_cog_tabs.append((fig, ax, canvas))

                    nested_tab_widget.addTab(sub_tab, f"{tab_names[j]}")

                # Store the plots specifically for X_cog
                self.plots[component] = x_cog_tabs

            else:
                # Initialize three plots for other components
                plots = []
                for j in range(3):
                    fig, ax = plt.subplots(figsize=(4, 3))  # Normal figure size
                    canvas = FigureCanvas(fig)
                    tab_layout.addWidget(canvas, 0, j)
                    plots.append((fig, ax, canvas))

                # Store the plots
                self.plots[component] = plots

            # Update plots based on the component
            self.update_plots(component)

            # Add the tab to the main tab widget
            self.tab_widget.addTab(tab, component)
            self.add_placeholder_icons(component)

    def update_plots(self, component):
        """Update the plots for the given component."""
        self.calculation_results = self.gui.calculation_results
        self.previous_calculation_results = self.gui.previous_calculation_results

        if component not in self.plots:
            print(f"Error: Component {component} not found in self.plots.")
            return

        if component not in self.calculation_results:
            print(f"Error: Component {component} not found in calculation results.")
            return

        for i, (fig, ax, canvas) in enumerate(self.plots[component]):
            ax.clear()
            if component in ["Duct", "Pylon", "Empennage", "Support", "Nacelle", "Control"]:
                self.update_aerodynamics_plot(component, i, fig, ax, canvas)
            elif component == "Inflow":
                self.update_inflow_plot(i, fig, ax, canvas)
            elif component == "Weight":
                self.update_weight_plot(i, fig, ax, canvas)
            elif component == "CD0":
                self.update_cd0_empennage(i, fig, ax, canvas)
            elif component in ["X_cog", "Vtail", "Htail"]:
                self.update_flight_mechanics_plot(component, i, fig, ax, canvas)

            elif component in ["Geometry"]:
                self.update_area_plot(i, fig, ax, canvas)

            canvas.draw()

        print(f"Updated plots for {component} tab")  # Debug print

    def update_aerodynamics_plot(self, component, plot_index, fig, ax, canvas):
        """Update the plot for Aerodynamics."""
        if component == 'Pylon':
            if plot_index == 0:
                ax.clear()
                file = 'pylon' + self.parameters['pylon_profile']
                alpha_pylon = self.calculation_results['Pylon']["Inflow"][0]
                polars_pylon = get_polar(file)

                # Left y-axis
                ax1 = ax
                ax1.plot(polars_pylon[0], polars_pylon[1], label=r'$C_l$', color='tab:blue')
                ax1.axvline(alpha_pylon, color='black', alpha=0.5, linestyle='dashed', label=r'$\alpha_{pylon}$')
                ax1.plot(alpha_pylon, self.calculation_results['Pylon']["Cl"][0], marker='o', color='tab:blue')

                ax1.set_ylabel(r'$C_l$ [-]', color='tab:blue')
                ax1.tick_params(axis='y', labelcolor='tab:blue')

                if not hasattr(ax, 'secondary_ax'):
                    ax.secondary_ax = ax.twinx()
                ax.secondary_ax.clear()

                ax.secondary_ax.plot(polars_pylon[0], polars_pylon[2], label=r'$C_d$', color='tab:red')
                ax.secondary_ax.plot(alpha_pylon, self.calculation_results['Pylon']["Cd"][0], marker='o', color='tab:red')
                ax.secondary_ax.plot(polars_pylon[0], polars_pylon[3], label=r'$C_m$', color='tab:green')
                ax.secondary_ax.plot(alpha_pylon, self.calculation_results['Pylon']["Cm"][0], marker='o', color='tab:green')
                ax.secondary_ax.set_ylabel(r'$C_d$, $C_m$ [-]', color='tab:red')
                ax.secondary_ax.yaxis.set_label_coords(1.20, 0.5)
                ax.secondary_ax.tick_params(axis='y', labelcolor='tab:red')

                # Common x-axis

                ax1.set_xlabel(r'$\alpha$ [deg]')
                ax1.set_title('XFoil Coefficients')
                ax1.grid(True)
                ax1.get_figure().tight_layout()

                # Combine legends from both axes
                lines_1, labels_1 = ax1.get_legend_handles_labels()
                lines_2, labels_2 =  ax.secondary_ax.get_legend_handles_labels()
                ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')
            elif plot_index == 1:
                ax.clear()
                ax.set_axis_off()
                ax.set_facecolor('white')  # Set the background color to white
                fig.patch.set_facecolor('white')
            elif plot_index == 2:
                columns = ['Parameter', 'Local value', 'Normalized value']
                cell_text = [
                    [r'$\alpha_{in}$ [deg]', f"{self.calculation_results['Pylon']['Inflow'][0]:.3f}", ' '],
                    [r'$V_{in}$ [m/s]', f"{self.calculation_results['Pylon']['Inflow'][1]:.3f}", ' '],
                    [r'$C_L$ [-]', f"{self.calculation_results['Pylon']['Cl'][0]:.6f}",
                     f"{self.calculation_results['Pylon']['Cl'][1]:.6f}"],
                    [r'$C_{D0}$ [-]', f"{0:.6f}", f"{0:.6f}"],
                    [r'$C_{Di}$ [-]', f"{0:.6f}", f"{0:.6f}"],
                    [r'$C_D$ [-]', f"{self.calculation_results['Pylon']['Cd'][0]:.6f}",
                     f"{self.calculation_results['Pylon']['Cd'][1]:.6f}"],
                    [r'$C_M$ [-]', f"{self.calculation_results['Pylon']['Cm'][0]:.6f}",
                     f"{self.calculation_results['Pylon']['Cm'][1]:.6f}"],
                    [r'Weight [kg]', f"{self.calculation_results['Weight']['w_vector_duuc'][3]:.2f}", ' ']
                ]

                style_table(ax, cell_text, columns, title='Pylon')

        elif component == 'Support':
            if plot_index == 0:
                ax.clear()
                file = 'support' + self.parameters['support_profile']
                alpha_support = self.calculation_results['Support']['Inflow'][0]
                polars_support = get_polar(file)

                # Left y-axis
                ax1 = ax
                ax1.plot(polars_support[0], polars_support[1], label=r'$C_l$', color='tab:blue')
                ax1.axvline(alpha_support, color='black', alpha=0.5, linestyle='dashed', label=r'$\alpha_{support}$')
                ax1.plot(alpha_support, self.calculation_results['Support']['Cl'][0], marker='o', color='tab:blue')

                ax1.set_ylabel(r'$C_l$ [-]', color='tab:blue')
                ax1.tick_params(axis='y', labelcolor='tab:blue')

                if not hasattr(ax, 'secondary_ax'):
                    ax.secondary_ax = ax.twinx()

                ax.secondary_ax.clear()
                ax.secondary_ax.plot(polars_support[0], polars_support[2], label=r'$C_d$', color='tab:red')
                ax.secondary_ax.plot(alpha_support, self.calculation_results['Support']['Cd'][0], marker='o',
                                     color='tab:red')
                ax.secondary_ax.plot(polars_support[0], polars_support[3], label=r'$C_m$', color='tab:green')
                ax.secondary_ax.plot(alpha_support, self.calculation_results['Support']['Cm'][0], marker='o',
                                     color='tab:green')
                ax.secondary_ax.set_ylabel(r'$C_d$, $C_m$ [-]', color='tab:red')
                ax.secondary_ax.yaxis.set_label_coords(1.20, 0.5)
                ax.secondary_ax.tick_params(axis='y', labelcolor='tab:red')

                # Common x-axis

                ax1.set_xlabel(r'$\alpha$ [deg]')
                ax1.set_title('XFoil Coefficients')
                ax1.grid(True)
                ax1.get_figure().tight_layout()

                # Combine legends from both axes
                lines_1, labels_1 = ax1.get_legend_handles_labels()
                lines_2, labels_2 = ax.secondary_ax.get_legend_handles_labels()
                ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')
            elif plot_index == 1:
                ax.clear()
                ax.set_axis_off()
                ax.set_facecolor('white')  # Set the background color to white
                fig.patch.set_facecolor('white')
            elif plot_index == 2:
                columns = ['Parameter', 'Local value', 'Normalized value']
                cell_text = [
                    [r'$\alpha_{in}$ [deg]', f"{self.calculation_results['Support']['Inflow'][0]:.3f}", ' '],
                    [r'$V_{in}$ [m/s]', f"{self.calculation_results['Support']['Inflow'][1]:.3f}", ' '],
                    [r'$C_L$ [-]', f"{self.calculation_results['Support']['Cl'][0]:.6f}",
                     f"{self.calculation_results['Support']['Cl'][1]:.6f}"],
                    [r'$C_{D0}$ [-]', f"{0:.6f}", f"{0:.6f}"],
                    [r'$C_{Di}$ [-]', f"{0:.6f}", f"{0:.6f}"],
                    [r'$C_D$ [-]', f"{self.calculation_results['Support']['Cd'][0]:.6f}",
                     f"{self.calculation_results['Support']['Cd'][1]:.6f}"],
                    [r'$C_M$ [-]', f"{self.calculation_results['Support']['Cm'][0]:.6f}",
                     f"{self.calculation_results['Support']['Cm'][1]:.6f}"],
                    [r'Weight [kg]', f"{self.calculation_results['Weight']['w_vector_duuc'][4]:.2f}", ' ']
                ]

                style_table(ax, cell_text, columns, title='Support')

        elif component == 'Duct':
            aspect_duct = self.parameters['duct_diameter'] / self.parameters['duct_chord']
            if plot_index == 0:
                ax.clear()
                alpha_duct = self.calculation_results['Duct']['Inflow'][0]
                a_vect = np.linspace(0, 15, 31)
                cl_vect = self.calculation_results["Empennage"]["Vectors"][5][:, 0]
                cd_vect = self.calculation_results["Empennage"]["Vectors"][5][:, 1]
                cm_vect = self.calculation_results["Empennage"]["Vectors"][5][:, 2]

                # Left y-axis
                ax1 = ax
                ax1.plot(a_vect, cl_vect, label=r'$C_l$', color='tab:blue')
                ax1.axvline(alpha_duct, color='black', alpha=0.5, linestyle='dashed', label=r'$\alpha_{duct}$')
                ax1.plot(alpha_duct, self.calculation_results['Duct']['Cl'][0], marker='o', color='tab:blue')

                ax1.set_ylabel(r'$C_l$ [-]', color='tab:blue')
                ax1.tick_params(axis='y', labelcolor='tab:blue')

                if not hasattr(ax, 'secondary_ax'):
                    ax.secondary_ax = ax.twinx()
                ax.secondary_ax.clear()
                ax.secondary_ax.plot(a_vect, cd_vect, label=r'$C_d$', color='tab:red')
                ax.secondary_ax.plot(alpha_duct, self.calculation_results['Duct']['Cd'][0], marker='o', color='tab:red')
                ax.secondary_ax.plot(a_vect, cm_vect, label=r'$C_m$', color='tab:green')
                ax.secondary_ax.plot(alpha_duct, self.calculation_results['Duct']['Cm'][0], marker='o', color='tab:green')
                ax.secondary_ax.set_ylabel(r'$C_d$, $C_m$ [-]', color='tab:red')
                ax.secondary_ax.yaxis.set_label_coords(1.15, 0.5)
                ax.secondary_ax.tick_params(axis='y', labelcolor='tab:red')

                # Common x-axi
                ax1.set_xlabel(r'$\alpha$ [deg]')
                ax1.set_title('Coefficients')
                ax1.grid(True)
                ax1.get_figure().tight_layout()

                # Combine legends from both axes
                lines_1, labels_1 = ax1.get_legend_handles_labels()
                lines_2, labels_2 = ax.secondary_ax.get_legend_handles_labels()
                ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')
            elif plot_index == 1:
                ax.clear()
                cl_vect = self.calculation_results["Empennage"]["Vectors"][5][:, 0]
                cd_vect = self.calculation_results["Empennage"]["Vectors"][5][:, 1]
                cd0_duct = self.calculation_results["CD0"]['cd0_duuc'][0]/2
                cl_the = []
                cd_the = []
                a_vect = np.linspace(0, 15, 31)

                for i in range(len(a_vect)):
                    al = np.radians(a_vect[i])
                    kp = 6.25 * np.sin(aspect_duct / 2)
                    kv = np.pi / 3
                    cl_theory = kp * np.sin(al) * np.cos(al) ** 2 + kv * np.cos(al) * np.sin(al) ** 2
                    cl_the.append(cl_theory)
                    cd_the.append(cd0_duct + 0.06 * cl_theory ** 2)
                # Left y-axis
                ax1 = ax
                ax1.plot(cd_vect, cl_vect, label=r'Duct', color='tab:blue')
                ax1.plot(cd_the, cl_the, label='Leading Edge Suction analogy', color="tab:orange", linestyle="--")
                ax1.plot([cd0_duct + cla ** 2 / (2 * np.pi * aspect_duct) for cla in ref5r.cla_cl_ar_1_52],
                         [cl * 0.95 for cl in ref5r.cla_cl_ar_1_52], label=r'Experimental (AR = 1.5)', color='tab:purple',
                         linestyle='dashed')
                ax1.plot(self.calculation_results['Duct']['Cd'][0], self.calculation_results['Duct']['Cl'][0], marker='o',
                         color='tab:blue')
                ax1.set_ylabel(r'$C_l$ [-]')
                ax1.set_xlabel(r'$C_d$ [-]')
                ax1.set_title(r'$C_l$ vs. $C_d$')
                ax1.grid(True)
                ax1.legend()
                ax1.get_figure().tight_layout()
            elif plot_index == 2:
                columns = ['Parameter', 'Local value', 'Normalized value']
                cell_text = [
                    [r'$\alpha_{in}$ [deg]', f"{self.calculation_results['Duct']['Inflow'][0]:.3f}", ' '],
                    [r'$V_{in}$ [m/s]', f"{self.calculation_results['Duct']['Inflow'][1]:.3f}", ' '],
                    [r'$C_L$ [-]', f"{self.calculation_results['Duct']['Cl'][0]:.6f}",
                     f"{self.calculation_results['Duct']['Cl'][1]:.6f}"],
                    [r'$C_D$ [-]', f"{self.calculation_results['Duct']['Cd'][0]:.6f}",
                     f"{self.calculation_results['Duct']['Cd'][1]:.6f}"],
                    [r'$C_M$ [-]', f"{self.calculation_results['Duct']['Cm'][0]:.6f}",
                     f"{self.calculation_results['Duct']['Cm'][1]:.6f}"],
                    [r'Weight [kg]', f"{self.calculation_results['Weight']['w_vector_duuc'][2]:.2f}", ' '],
                    [r'Aspect Ratio [-]', aspect_duct, '']
                ]
                style_table(ax, cell_text, columns, title='Duct')

                reference_text = (r"$\bf{Reference\ Data\ Source:}$" "\n"
                                  r"Leading Edge Suction Analogy - Maqsood 2013" "\n\n"
                                  r"Experimental Investigation of Lift, Drag, and Pitching Moment" 
                                  "\nof Five Annular Airfoils - Fletcher 1957")

                ax.text(0, .12, reference_text, horizontalalignment='left',
                        verticalalignment='center', transform=ax.transAxes, fontsize=12,
                        bbox=dict(facecolor='lightblue', edgecolor='black', alpha=0.7))

        elif component == 'Control':
            if plot_index == 0:
                ax.clear()
                file = 'control' + self.parameters['cv_airfoil']
                alpha_elev = self.calculation_results['Control']['Inflow'][0]
                alpha_rud = self.calculation_results['Control']['Inflow'][2]
                polars_elev = get_polar(file)

                # Left y-axis
                ax1 = ax
                ax1.plot(polars_elev[0], polars_elev[1], label=r'$C_l$', color='tab:blue')
                ax1.axvline(alpha_elev, color='black', alpha=0.5, linestyle='dashed', label=r'$\alpha_{elevator}$')
                ax1.axvline(alpha_rud, color='black', alpha=0.5, linestyle='dashed', label=r'$\alpha_{rudder}$')
                ax1.plot(alpha_elev, self.calculation_results['Control']['Cl'][0], marker='o', color='tab:blue')
                ax1.plot(alpha_rud, self.calculation_results['Control']['Cy'][0], marker='o', color='tab:purple',
                         label='Rudder')

                ax1.set_ylabel(r'$C_l$ [-]', color='tab:blue')
                ax1.tick_params(axis='y', labelcolor='tab:blue')

                if not hasattr(ax, 'secondary_ax'):
                    ax.secondary_ax = ax.twinx()
                ax.secondary_ax.clear()
                ax.secondary_ax.plot(polars_elev[0], polars_elev[2], label=r'$C_d$', color='tab:red')
                ax.secondary_ax.plot(alpha_elev, self.calculation_results['Control']['Cd'][0], marker='o', color='tab:red')
                ax.secondary_ax.plot(alpha_elev, self.calculation_results['Control']['Cd'][2], marker='o',
                                     color='tab:orange', label='Rudder')
                ax.secondary_ax.plot(polars_elev[0], polars_elev[3], label=r'$C_m$', color='tab:green')
                ax.secondary_ax.plot(alpha_elev, self.calculation_results['Control']['Cm'][0], marker='o', color='tab:green')
                ax.secondary_ax.plot(alpha_elev, self.calculation_results['Control']['Cm'][2], marker='o',
                                     color='tab:grey', label='Rudder')
                ax.secondary_ax.set_ylabel(r'$C_d$, $C_m$ [-]', color='tab:red')
                ax.secondary_ax.yaxis.set_label_coords(1.20, 0.5)
                ax.secondary_ax.tick_params(axis='y', labelcolor='tab:red')

                # Common x-axis
                ax1.set_xlabel(r'$\alpha$ [deg]')
                ax1.set_title('XFoil Coefficients')
                ax1.grid(True)
                ax1.get_figure().tight_layout()

                # Combine legends from both axes
                lines_1, labels_1 = ax1.get_legend_handles_labels()
                lines_2, labels_2 = ax.secondary_ax.get_legend_handles_labels()
                ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')
            elif plot_index == 1:
                ax.clear()
                ax.set_axis_off()
                ax.set_facecolor('white')  # Set the background color to white
                fig.patch.set_facecolor('white')
            elif plot_index == 2:
                columns = ['Parameter', 'Local value', 'Normalized value']
                cell_text = [
                    [r'$\alpha_{in-elevator}$ [deg]', f"{self.calculation_results['Control']['Inflow'][0]:.3f}", ' '],
                    [r'$V_{in-elevator}$ [m/s]', f"{self.calculation_results['Control']['Inflow'][1]:.3f}", ' '],
                    [r'$\alpha_{in-rudder}$ [deg]', f"{self.calculation_results['Control']['Inflow'][2]:.3f}", ' '],
                    [r'$V_{in-rudder}$ [m/s]', f"{self.calculation_results['Control']['Inflow'][3]:.3f}", ' '],
                    [r'$C_{L-elevator}$ [-]', f"{self.calculation_results['Control']['Cl'][0]:.6f}",
                     f"{self.calculation_results['Control']['Cl'][1]:.6f}"],
                    [r'$C_{Y-rudder}$ [-]', f"{self.calculation_results['Control']['Cy'][0]:.6f}",
                     f"{self.calculation_results['Control']['Cy'][0]:.6f}"],
                    [r'$C_{D-elevator}$ [-]', f"{self.calculation_results['Control']['Cd'][0]:.6f}",
                     f"{self.calculation_results['Control']['Cd'][1]:.6f}"],
                    [r'$C_{D-rudder}$ [-]', f"{self.calculation_results['Control']['Cd'][2]:.6f}",
                     f"{self.calculation_results['Control']['Cd'][3]:.6f}"],
                    [r'$C_{M-elevator}$ [-]', f"{self.calculation_results['Control']['Cm'][0]:.6f}",
                     f"{self.calculation_results['Control']['Cm'][1]:.6f}"],
                    [r'$C_{M-rudder}$ [-]', f"{self.calculation_results['Control']['Cm'][2]:.6f}",
                     f"{self.calculation_results['Control']['Cm'][3]:.6f}"],
                    [r'Weight [kg]', f"{self.calculation_results['Weight']['w_vector_duuc'][5]:.2f}", ' ']
                ]

                style_table(ax, cell_text, columns, title='Rudder')

        elif component == 'Nacelle':
            if plot_index == 0:
                ax.clear()
                alpha_nacelle = self.calculation_results['Nacelle']['Inflow'][0]
                a_vector = np.linspace(-5, 15, 21)

                # Left y-axis
                ax1 = ax
                ax1.plot(a_vector, self.calculation_results['Nacelle']['Vector'][0], label=r'$C_l$', color='tab:blue')
                ax1.axvline(alpha_nacelle, color='black', alpha=0.5, linestyle='dashed', label=r'$\alpha_{nacelle}$')
                ax1.plot(alpha_nacelle, self.calculation_results['Nacelle']['Cl'][0], marker='o', color='tab:blue')
                ax1.plot(a_vector, self.calculation_results['Nacelle']['Vector'][2], label=r'$C_m$', color='tab:green')
                ax1.plot(alpha_nacelle, self.calculation_results['Nacelle']['Cm'][0], marker='o', color='tab:green')

                ax1.set_ylabel(r'$C_l$, $C_m$ [-]', color='tab:blue')
                ax1.tick_params(axis='y', labelcolor='tab:blue')

                # Right y-axis (secondary axis)
                if not hasattr(ax, 'secondary_ax'):
                    ax.secondary_ax = ax.twinx()

                # Clear the secondary axis explicitly
                ax.secondary_ax.clear()

                ax.secondary_ax.plot(a_vector, self.calculation_results['Nacelle']['Vector'][1], label=r'$C_d$', color='tab:red')
                ax.secondary_ax.plot(alpha_nacelle, self.calculation_results['Nacelle']['Cd'][0], marker='o', color='tab:red')

                ax.secondary_ax.set_ylabel(r'$C_d$ [-]', color='tab:red')
                ax.secondary_ax.yaxis.set_label_coords(1.20, 0.5)
                ax.secondary_ax.tick_params(axis='y', labelcolor='tab:red')

                # Common x-axis
                ax1.set_xlabel(r'$\alpha$ [deg]')
                ax1.set_title('Coefficients')
                ax1.grid(True)
                ax1.get_figure().tight_layout()

                # Combine legends from both axes
                lines_1, labels_1 = ax1.get_legend_handles_labels()
                lines_2, labels_2 = ax.secondary_ax.get_legend_handles_labels()
                ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')
            elif plot_index == 1:
                ax.clear()
                ax.set_axis_off()
                ax.set_facecolor('white')  # Set the background color to white
                fig.patch.set_facecolor('white')
            elif plot_index == 2:
                columns = ['Parameter', 'Local value', 'Normalized value']
                cell_text = [
                    [r'$\alpha_{in}$ [deg]', f"{self.calculation_results['Nacelle']['Inflow'][0]:.3f}", ' '],
                    [r'$V_{in}$ [m/s]', f"{self.calculation_results['Nacelle']['Inflow'][1]:.3f}", ' '],
                    [r'$C_L$ [-]', f"{self.calculation_results['Nacelle']['Cl'][0]:.6f}",
                     f"{self.calculation_results['Nacelle']['Cl'][1]:.6f}"],
                    [r'$C_D$ [-]', f"{self.calculation_results['Nacelle']['Cd'][0]:.6f}",
                     f"{self.calculation_results['Nacelle']['Cd'][1]:.6f}"],
                    [r'$C_M$ [-]', f"{self.calculation_results['Nacelle']['Cm'][0]:.6f}",
                     f"{self.calculation_results['Nacelle']['Cm'][1]:.6f}"],
                    [r'Weight [kg]', f"{self.calculation_results['Weight']['w_vector_duuc'][7]:.2f}", ' ']
                ]

                style_table(ax, cell_text, columns, title='Nacelle')

        elif component == 'Empennage':
            a_vector = np.linspace(0, 15, 31)
            if plot_index == 0:
                ax.clear()
                cl_vector = self.calculation_results["Empennage"]["Vectors"][0]
                cd_vector = self.calculation_results["Empennage"]["Vectors"][1]
                cm_vector = self.calculation_results["Empennage"]["Vectors"][2]
                cl = self.calculation_results["Empennage"]["Cl"][0]
                cd = self.calculation_results["Empennage"]["Cd"][0]
                cm = self.calculation_results["Empennage"]["Cm"][0]
                alpha = self.calculation_results["Empennage"]["Inflow"][0]

                # Left y-axis
                ax1 = ax
                ax1.plot(a_vector, cl_vector, label=r'$C_l$', color='tab:blue')
                ax1.axvline(alpha, color='black', alpha=0.5, linestyle='dashed', label=r'$\alpha_{empennage}$')
                ax1.plot(alpha, cl, marker='o', color='tab:blue')

                ax1.set_ylabel(r'$C_l$ [-]', color='tab:blue')
                ax1.tick_params(axis='y', labelcolor='tab:blue')

                ax1.plot(a_vector, cd_vector, label=r'$C_d$', color='tab:red')
                ax1.plot(alpha, cd, marker='o', color='tab:red')
                ax1.plot(a_vector, cm_vector, label=r'$C_m$', color='tab:green')
                ax1.plot(alpha, cm, marker='o', color='tab:green')
                ax1.tick_params(axis='y', labelcolor='tab:red')

                # Common x-axis
                ax1.set_xlabel(r'$\alpha$ [deg]')
                ax1.set_title('Coefficients 1 Empennage')
                ax1.grid(True)
                ax1.get_figure().tight_layout()

                ax1.legend(loc='best')
            elif plot_index == 1:
                cl_duct = self.calculation_results["Empennage"]["Vectors"][5][:, 3]
                cl_pylon = self.calculation_results["Empennage"]["Vectors"][6][:, 3]
                cl_support = self.calculation_results["Empennage"]["Vectors"][7][:, 3]
                cl_control = self.calculation_results["Empennage"]["Vectors"][8][:, 3]
                cl_nacelle = self.calculation_results["Empennage"]["Vectors"][9][:, 3]

                ax.clear()
                ax.plot(a_vector, cl_duct, label=r"Duct")
                ax.plot(a_vector, cl_pylon, label="Pylon")
                ax.plot(a_vector, cl_support, label="Support")
                ax.plot(a_vector, cl_control, label="Control")
                ax.plot(a_vector, cl_nacelle, label="Nacelle")
                ax.axvline(self.calculation_results["Empennage"]["Inflow"][0], color='black', alpha=0.5,
                           linestyle='dashed')
                ax.set_title("Lift Breakdown Propulsive Empennage")
                ax.set_ylabel(r"$C_L$ [-]")
                ax.set_xlabel(r"$\alpha$ [deg]")
                ax.grid(True)
                ax.legend()
            elif plot_index == 2:
                cd_duct = self.calculation_results["Empennage"]["Vectors"][5][:, 4]
                cd_pylon = self.calculation_results["Empennage"]["Vectors"][6][:, 4]
                cd_support = self.calculation_results["Empennage"]["Vectors"][7][:, 4]
                cd_control = self.calculation_results["Empennage"]["Vectors"][8][:, 4]
                cd_nacelle = self.calculation_results["Empennage"]["Vectors"][9][:, 4]

                ax.clear()
                ax.plot(a_vector, cd_duct, label=r"Duct")
                ax.plot(a_vector, cd_pylon, label="Pylon")
                ax.plot(a_vector, cd_support, label="Support")
                ax.plot(a_vector, cd_control, label="Control")
                ax.plot(a_vector, cd_nacelle, label="Nacelle")
                ax.axvline(self.calculation_results["Empennage"]["Inflow"][0], color='black', alpha=0.5,
                           linestyle='dashed')
                ax.set_title("Drag Breakdown Propulsive Empennage")
                ax.set_ylabel(r"$C_D$ [-]")
                ax.set_xlabel(r"$\alpha$ [deg]")
                ax.grid(True)
                ax.legend()

    def update_flight_mechanics_plot(self, component, plot_index, fig, ax, canvas):
        """Update the plot for Flight Mechanics."""
        if component == "Weight":
            w_vector_duuc = self.calculation_results['Weight']['w_vector_duuc']
            w_vector_atr72_600 = self.calculation_results['Weight']['w_vector_atr']

            prev_w_vector_duuc = self.previous_calculation_results['Weight'][
                'w_vector_duuc'] if self.previous_calculation_results else None
            prev_w_vector_atr72_600 = self.previous_calculation_results['Weight'][
                'w_vector_atr'] if self.previous_calculation_results else None

            self.plot_weight_distribution(w_vector_duuc, w_vector_atr72_600, prev_w_vector_duuc,
                                          prev_w_vector_atr72_600, plot_index, ax)
        elif component == "X_cog":
            self.update_x_cog_plot(plot_index, fig, ax, canvas)
            ax.set_title("Center of Gravity Positions")

        elif component == "Vtail":
            current_data = self.calculation_results['Vtail']
            prev_data = self.previous_calculation_results['Vtail'] if self.previous_calculation_results else None

            if plot_index == 0:
                self.plot_cy_vs_sv(ax, current_data, prev_data)
            elif plot_index == 1:
                self.plot_cybeta_vs_sv(ax, current_data, prev_data)
            elif plot_index == 2:
                # Empty plot or additional information
                ax.set_axis_off()
                ax.set_facecolor('white')
                fig.patch.set_facecolor('white')

        elif component == "Htail":
            self.plot_xplot(ax, plot_index)

        canvas.draw()

    def plot_xplot(self, ax, plot_index):
        a1_atr, b1_atr, a2_atr = self.calculation_results["Htail"]["ATR"]
        a1_duuc, b1_duuc, a2_duuc, stm = self.calculation_results["Htail"]["DUUC"]
        x_cg_duuc = (self.calculation_results["X_cog"]["x_cog_duuc"][2]
                     - self.calculation_results["X_cog"]["x_cog_duuc"][3] - 0.25)
        x_cg_atr = (self.calculation_results["X_cog"]["x_cog_atr"][2]
                    - self.calculation_results["X_cog"]["x_cog_atr"][3] - 0.25)

        cmac = 2.2345  # Mean aerodynamic chord
        s_w = 61  # Wing area (update this based on real reference if needed)
        static_margin = stm / 100  # Convert to fraction

        def compute_xplot_data(a1, b1, a2, xcg, color, label):
            x1 = np.linspace(-1, 0, 101)
            x2 = np.linspace(-1, 2, 101)
            x_ac = 0.25 * cmac

            x_cg_ac = x1 + x_ac
            y1 = a1 * x_cg_ac * cmac + b1
            y2 = a2 * x2
            y3 = y2 + static_margin

            y_intersection = a2 * xcg + static_margin
            x_intersect = (y_intersection - b1) / (a1 * cmac) - x_ac

            sh_sw = np.round(y_intersection, 4)
            new_sh = y_intersection * s_w
            cg_range = xcg - x_intersect

            # Plotting
            ax.plot(x1, y1, label=f'{label} Control', color=color)
            ax.plot(x2, y2, color=color, linestyle='-', alpha=0.4)
            ax.plot(x2, y3, label=f'{label} Static Margin', color=color, linestyle='--')
            ax.vlines(xcg, 0, y_intersection, color=color, linestyle="dotted")
            ax.scatter([xcg, x_intersect], [y_intersection, y_intersection], color=color, zorder=3)
            ax.hlines(y_intersection, x_intersect, xcg, color=color, linestyle="dotted")
            return new_sh

        if plot_index == 0:
            ax.clear()

            compute_xplot_data(a1_atr, b1_atr, a2_atr, x_cg_atr, color="tab:orange", label="ATR")
            compute_xplot_data(a1_duuc, b1_duuc, a2_duuc, x_cg_duuc, color="tab:blue", label="DUUC")

            # Set plot labels and styling
            ax.set_xlabel('Center of Gravity [m]')
            ax.set_ylabel(r'$S_H/S_W$ [-]')
            ax.set_title("X-Plot: ATR vs DUUC")
            ax.set_ylim([-0.6, 0.8])
            ax.grid(True)
            ax.legend()

    def plot_cy_vs_sv(self, ax, current_data, prev_data):
        s_array, s_array_atr, s_stab_duuc, s_stab_atr, a, b, cyd = current_data
        array = np.linspace(0.1, 3.0, len(s_array))  # Assuming 301 points as in your original setup

        ax.plot(array, s_array, label=r'Prediction line DUUC - control', color="tab:blue")
        ax.plot(array, s_array_atr, label=r'Prediction line ATR - control', color="tab:orange")

        # ATR value
        ax.plot([0, 0.975], [a, a], color="tab:orange", linestyle="dashed")
        ax.plot(0.975, a, 'o', markersize=6, color="tab:orange", label="ATR - value")
        ax.plot([0.975, 0.975], [0, a], linestyle='dashed', color="tab:orange")

        # Current DUUC value
        ax.plot([0, cyd], [b, b], color="tab:blue", linestyle="dashed")
        ax.plot(cyd, b, 'o', markersize=6, color="tab:blue", label="DUUC - current")
        ax.plot([cyd, cyd], [0, b], color="tab:blue", linestyle="dashed")

        # Previous DUUC value (if available)
        if prev_data:
            prev_s_array, prev_s_array_atr, _, _, prev_a, prev_b, prev_cyd = prev_data
            ax.plot([0, prev_cyd], [prev_b, prev_b], color="tab:blue", linestyle="dashed", alpha=0.3)
            ax.plot(prev_cyd, prev_b, 'x', markersize=6, color="tab:blue", alpha=0.3, label="DUUC - previous")
            ax.plot([prev_cyd, prev_cyd], [0, prev_b], color="tab:blue", linestyle="dashed", alpha=0.3)

        ax.set_xlim(0, 2)
        ax.set_ylim(0, 50)
        ax.set_xlabel(r'$C_{Y}$ [-]')
        ax.set_ylabel(r'$S_{V}$ [$m2$]')
        ax.set_title(r'$S_{V}$ vs. $C_{Y}$')
        ax.legend()
        ax.grid(True)

    def plot_cybeta_vs_sv(self, ax, current_data, prev_data):
        s_array, s_array_atr, s_stab_duuc, s_stab_atr, a, b, cyd = current_data
        array = np.linspace(0.1, 3.0, len(s_stab_duuc))

        # Plot previous iteration's line with reduced opacity
        if prev_data:
            _, _, prev_s_stab_duuc, _, _, _, _ = prev_data
            ax.plot(array, prev_s_stab_duuc, label='Previous DUUC - stability', color="tab:blue", alpha=0.3)

        # Plot current iteration's line with full opacity
        ax.plot(array, s_stab_duuc, label='Current DUUC - stability', color="tab:blue")
        ax.plot(array, s_stab_atr, label='ATR - stability', color="tab:orange")

        # ATR value
        ax.plot([0, 2.228], [15.01, 15.01], color="tab:orange", linestyle="dashed")
        ax.plot(2.228, 15.01, 'o', markersize=6, color="tab:orange", label="ATR - value")
        ax.plot([2.228, 2.228], [0, 15.01], linestyle='dashed', color="tab:orange")

        ax.set_xlim(1, 3)
        ax.set_ylim(0, 50)
        ax.set_xlabel(r'$C_{Y_{\beta}}$ [-]')
        ax.set_ylabel(r'$S_{V}$ [$m2$]')
        ax.set_title(r'$S_{V}$ vs. $C_{Y_{\beta}}$')
        ax.legend()
        ax.grid(True)

    def update_weight_plot(self, plot_index, fig, ax, canvas):
        w_vector_duuc = self.calculation_results['Weight']['w_vector_duuc']
        w_vector_atr = self.calculation_results['Weight']['w_vector_atr']

        prev_w_vector_duuc = self.previous_calculation_results['Weight'][
            'w_vector_duuc'] if self.previous_calculation_results else None
        prev_w_vector_atr = self.previous_calculation_results['Weight'][
            'w_vector_atr'] if self.previous_calculation_results else None

        self.plot_weight_distribution(w_vector_duuc, w_vector_atr, prev_w_vector_duuc, prev_w_vector_atr,
                                      plot_index, ax)

    def update_inflow_plot(self, plot_index, fig, ax, canvas):
        """Update the inflow plot."""
        v_input = self.calculation_results['Inflow']['v_inflow']
        a_input = self.calculation_results['Inflow']['a_inflow']
        station = self.calculation_results['Inflow']['station']
        pylon = self.calculation_results['Inflow']['pylon_in']
        alpha = self.parameters['alpha']
        v_inf = self.parameters['velocity']
        d_duct = self.parameters['duct_diameter']
        naca = self.parameters['duct_profile']
        chord = self.parameters['duct_chord']
        d_prop = self.parameters['propeller_diameter']
        x_cv = self.parameters['x_control_vanes']
        x_sup = self.parameters['x_support']
        d_sup = area_ratio(naca, chord, d_duct/2, x_sup)[0]
        d_cv = area_ratio(naca, chord, d_duct/2, x_cv)[0]

        va_in = self.parameters["velocity"] * (0.25 * np.pi * self.parameters["duct_diameter"] ** 2)

        v_prop = va_in / (0.25 * np.pi * d_prop ** 2)
        v_sup = va_in / (0.25 * np.pi * d_sup ** 2)
        v_cv = va_in / (0.25 * np.pi * d_cv ** 2)

        v_input2 = [v_inf, v_inf, v_prop, v_prop, v_sup, v_cv, (v_cv + v_inf) / 2]
        a_input2 = [alpha, alpha, alpha / 2, 0, 0, 0, alpha/2]

        # Load background image
        img = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\propulsive_empennage.png")

        if plot_index == 0:
            # Plot for inflow velocity
            ax.clear()
            v_in = np.array(v_input)
            up_lim = v_in.max()
            low_lim = v_in.min()
            ax.imshow(img, extent=[0.35, 6.00, low_lim * 0.75, up_lim * 1.25], aspect='auto')
            ax.plot(station, v_input, label=r"$V_{inflow}$ - power on", color="tab:blue", marker='o')
            ax.plot(station, v_input2, label=r"$V_{inflow}$ - power off", color="tab:blue", linestyle="dotted",
                    marker='o')
            ax.plot([0, 2.55], [v_inf, pylon[0]], label=r'$V_{inflow} - pylon$', color='tab:blue', linestyle='dashed', marker='x')
            ax.set_xlim(0, 6)
            ax.set_ylim(low_lim * 0.75, up_lim * 1.25)
            ax.xaxis.set_major_locator(ticker.NullLocator())
            ax.set_ylabel(r'$V_{inflow}$ [m/s]')
            ax.set_title(r'Inflow Velocity per Component')
            ax.tick_params(axis='y')
            ax.grid(True)
            ax.legend(loc="upper left")

        elif plot_index == 1:
            # Plot for inflow angle
            ax.clear()
            a_in = np.array(a_input)
            up_lim = a_in.max()
            low_lim = a_in.min()
            ax.imshow(img, extent=[0.35, 6.00, -low_lim * 0.75, up_lim * 1.25], aspect='auto')
            ax.plot(station, a_input, label=r"$\alpha$ - power on", color="tab:orange", marker='o')
            ax.plot(station, a_input2, label=r"$\alpha$ - power off", color="tab:orange", linestyle="dotted",
                    marker='o')
            ax.plot([0, 2.55], [alpha, pylon[1]], label=r'$\alpha - pylon$', color='tab:orange', linestyle='dashed',
                    marker='x')
            ax.set_xlim(0, 6)
            ax.set_ylim(low_lim * 1.25, up_lim * 1.25)
            ax.xaxis.set_major_locator(ticker.NullLocator())
            ax.set_ylabel(r"$\alpha_{inflow}$ [deg]")
            ax.set_title(r'Inflow Angle per Component')
            ax.tick_params(axis='y')
            ax.grid(True)
            ax.legend(loc="upper left")

        elif plot_index == 2:
            # Empty plot for the third tab
            ax.clear()
            ax.set_axis_off()
            ax.set_facecolor('white')  # Set the background color to white
            fig.patch.set_facecolor('white')

        canvas.draw()

    def update_area_plot(self, plot_index, fig, ax, canvas):
        profile_duct = self.calculation_results['Geometry'][0]
        profile_pylon = self.calculation_results['Geometry'][1]
        profile_support = self.calculation_results['Geometry'][2]
        profile_control = self.calculation_results['Geometry'][3]

        coordinates_folder = r"C:\Users\tomva\pythonProject\DUUC\data\airfoil_coordinates"
        duct_file = os.path.join(coordinates_folder, "Naca" + profile_duct + ".txt")
        pylon_file = os.path.join(coordinates_folder, "Naca" + profile_pylon + ".txt")
        support_file = os.path.join(coordinates_folder, "Naca" + profile_support + ".txt")
        control_file = os.path.join(coordinates_folder, "Naca" + profile_control + ".txt")

        duct_data = np.loadtxt(duct_file)
        pylon_data = np.loadtxt(pylon_file)
        support_data = np.loadtxt(support_file)
        control_data = np.loadtxt(control_file)

        duct_x, duct_y = duct_data[:, 0], duct_data[:, 1]
        pylon_x, pylon_y = pylon_data[:, 0], pylon_data[:, 1]
        support_x, support_y = support_data[:, 0], support_data[:, 1]
        control_x, control_y = control_data[:, 0], control_data[:, 1]

        """Update the Area vs Chord plot."""
        if plot_index == 0:
            # Data preparation for Area vs. Chord
            x_c = np.linspace(0, 1, 101)
            area = []
            area_av = []
            duct_radius = self.parameters['duct_diameter'] / 2

            for i in range(len(x_c)):
                area.append(np.pi * area_ratio(self.parameters['duct_profile'], 1, duct_radius, x_c[i])[1] ** 2)
                area_av.append(
                    np.pi * area_ratio(self.parameters['duct_profile'], 1, duct_radius, x_c[i])[1] ** 2 -
                    cross_sectional_area(x_c[i], self.parameters['nacelle_diameter'], self.parameters['x_prop'])[0])

            # Clear the axes
            ax.clear()

            # Plot Area vs. Chord
            ax.plot(x_c, area, label=r'Duct area', color="tab:blue")
            ax.plot(x_c, area_av, label=r'Duct minus Nacelle area', color="tab:orange")
            ax.axvline(x=self.parameters["x_prop"], linestyle='--', color='black', label='Propeller location')

            # Add labels and title
            ax.set_xlabel(r'x/c [-]')
            ax.set_ylabel(r'Cross sectional area [$m^2$]')
            ax.set_title(r'Area vs. chord')

            # Add legend and grid
            ax.legend()
            ax.grid(True)

        if plot_index == 1:
            # Data preparation for Radius Duct vs Radius Available
            x_c = np.linspace(0, 1, 101)
            r_r1 = []
            r_r2 = []
            duct_radius = self.parameters['duct_diameter'] / 2

            for i in range(len(x_c)):
                r_r1.append(0.5 * area_ratio(self.parameters['duct_profile'], 1, duct_radius, x_c[i])[0] / duct_radius)
                r_r2.append((0.5 * area_ratio(self.parameters['duct_profile'], 1, duct_radius, x_c[i])[0] -
                             cross_sectional_area(x_c[i], self.parameters['nacelle_diameter'],
                                                  self.parameters["x_prop"])[1]) /
                            duct_radius)

            # Clear the axes
            ax.clear()

            # Plot Radius Duct vs Radius Available
            ax.plot(x_c, r_r1, label=r'Duct radius', color="tab:blue")
            ax.plot(x_c, r_r2, label=r'Duct minus Nacelle radius', color="tab:orange")

            # Add horizontal and vertical lines
            ax.axvline(x=self.parameters["x_prop"], linestyle='--', color='black', alpha=0.5, label='Propeller location')
            ax.axhline(y=1.0, linestyle='--', color='tab:blue', alpha=0.5)

            # Add labels and title
            ax.set_xlabel(r'x/c [-]')
            ax.set_ylabel(r'$R_{i}$ / $R_{center}$ [-]')
            ax.set_title(r'Radius duct vs Radius available')

            # Add legend and grid
            ax.legend()
            ax.grid(True)

        elif plot_index == 2:
            ax.clear()
            ax.plot(duct_x, duct_y, label=f"Duct Airfoil (NACA{profile_duct})")
            ax.plot(pylon_x, pylon_y, linestyle="dashed", label=f"Pylon Airfoil (NACA{profile_pylon})")
            ax.plot(support_x, support_y, linestyle="dotted", label=f"Support Airfoil (NACA{profile_support})")
            ax.plot(control_x, control_y, linestyle="dashed", label=f"Control Airfoil (NACA{profile_control})")
            ax.legend()
            ax.set_title("Airfoils Empennage")
            ax.set_ylabel("y/c [-]")
            ax.set_ylim([-0.20, 0.20])
            ax.set_xlabel("x/c [-]")
            ax.grid(True)
            plt.tight_layout()

        canvas.draw()

    def update_x_cog_plot(self, plot_index, fig, ax, canvas):
        """Update X_cog plot with horizontal position lines for both DUUC and ATR."""
        from data.read_data import read_text_file
        current_duuc = self.calculation_results['X_cog']['x_cog_duuc']
        current_atr = self.calculation_results['X_cog']['x_cog_atr']
        prev_duuc = prev_atr = None

        if self.previous_calculation_results and 'X_cog' in self.previous_calculation_results:
            prev_duuc = self.previous_calculation_results['X_cog']['x_cog_duuc']
            prev_atr = self.previous_calculation_results['X_cog']['x_cog_atr']

        categories = ['Fuselage CG', 'Wing CG', 'Overall CG', 'LEMAC']
        y_pos = np.arange(len(categories))
        offset = 0.2  # Vertical offset between DUUC and ATR lines

        ax.clear()

        if plot_index == 0:

            for i, (duuc_val, atr_val) in enumerate(zip(current_duuc, current_atr)):
                ax.plot([0, duuc_val], [y_pos[i] + offset, y_pos[i] + offset],
                        color='tab:blue', linewidth=2, label="DUUC")
                ax.plot(duuc_val, y_pos[i] + offset, marker='o', markersize=8, color='tab:blue')
                ax.text(duuc_val + 0.1, y_pos[i] + offset, f"{duuc_val:.2f}", va='center', ha='left', fontsize=10,
                        color='tab:blue')

                ax.plot([0, atr_val], [y_pos[i] - offset, y_pos[i] - offset],
                        color='tab:orange', linewidth=2, label="ATR")
                ax.plot(atr_val, y_pos[i] - offset, marker='o', markersize=8, color='tab:orange')
                ax.text(atr_val + 0.1, y_pos[i] - offset, f"{atr_val:.2f}", va='center', ha='left', fontsize=10,
                        color='tab:orange')

                if prev_duuc and prev_atr:
                    ax.plot([0, prev_duuc[i]], [y_pos[i] + offset - 0.05,
                                                y_pos[i] + offset - 0.05],
                            color='tab:blue', alpha=0.3, linewidth=2)

                    ax.plot([0, prev_atr[i]], [y_pos[i] - offset + 0.05,
                                               y_pos[i] - offset + 0.05],
                            color='tab:orange', alpha=0.3, linewidth=2)

            ax.set_yticks(y_pos)
            ax.set_yticklabels(categories)
            ax.set_xlabel('X Position [m]')
            ax.grid(True, axis='x')
            ax.legend(loc='upper left')
            ax.set_title('Center of Gravity Positions')

        elif plot_index == 1:
            fig.clf()
            ax = fig.add_subplot(111)
            ATR = read_text_file(r'C:\Users\tomva\pythonProject\DUUC\data\CG_breakdown\ATR_wing_group_cg_breakdown.txt')
            DUUC = read_text_file(
                r'C:\Users\tomva\pythonProject\DUUC\data\CG_breakdown\DUUC_wing_group_cg_breakdown.txt')
            ax.text(-0.1, 0.5, ATR, ha='left', va='center', fontsize=8)
            ax.axis('off')
            ax.text(0.60, 0.5, DUUC, ha='left', va='center', fontsize=8)
            ax.axis('off')
            ax.set_title("Wing group breakdown")
            canvas.draw()

        elif plot_index == 2:
            fig.clf()
            ax = fig.add_subplot(111)
            ATR = read_text_file(r'C:\Users\tomva\pythonProject\DUUC\data\CG_breakdown\ATR_fuselage_group_cg_breakdown.txt')
            DUUC = read_text_file(
                r'C:\Users\tomva\pythonProject\DUUC\data\CG_breakdown\DUUC_fuselage_group_cg_breakdown.txt')
            ax.text(-0.1, 0.5, ATR, ha='left', va='center', fontsize=8)
            ax.axis('off')
            ax.text(0.60, 0.5, DUUC, ha='left', va='center', fontsize=8)
            ax.axis('off')
            ax.set_title("Fuselage group breakdown")
            canvas.draw()
        elif plot_index == 3:
            l_fuselage = self.parameters["fuselage_length"]
            x_cg = self.calculation_results["requirements"]["x_cg"][0]
            x_ac_wing = self.parameters["x_wing"] + 0.25 * self.parameters['wing_c_root']
            x = np.linspace(0, l_fuselage, 100)

            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])  # Left is 3x wider than right

            # Left subplot
            fig.clear()  # clear the full figure
            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
            ax0 = fig.add_subplot(gs[0])
            ax1 = fig.add_subplot(gs[1])
            image = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\DUUC_side_view_tailless.png")

            # Plot background image (aircraft side view)
            ax0.imshow(image, extent=[0, l_fuselage, 0, l_fuselage], aspect='auto')

            ax0.plot(x_cg, x, label=r"$x_{cg-DUUC}$", color="tab:blue")
            ax0.axvline(min(x_cg), color="tab:blue", linestyle="dashed")
            ax0.axvline(max(x_cg), color="tab:blue", linestyle="dashed")
            ax0.axvspan(min(x_cg), max(x_cg), color="tab:blue", alpha=0.2)

            ax0.axvline(11.5586, color="tab:orange", linestyle="dashed", label="$x_{cg-atr}$")
            ax0.axvline(x_ac_wing, color="tab:red", label=r"$x_{ac-wing}$")

            ax0.plot([min(x_cg), l_fuselage], [0, 0], color="black", linestyle="dashed")
            ax0.plot(min(x_cg), 0, color="tab:blue", marker="o")
            ax0.plot([np.mean(x_cg), l_fuselage], [0.5 * l_fuselage, 0.5 * l_fuselage], color="black",
                     linestyle="dashed")
            ax0.plot(np.mean(x_cg), 0.5 * l_fuselage, color="tab:blue", marker="o")
            ax0.plot([max(x_cg), l_fuselage], [l_fuselage, l_fuselage], color="black", linestyle="dashed")
            ax0.plot(max(x_cg), l_fuselage, color="tab:blue", marker="o")

            ax0.set_title("Center of Gravity Shift with Position of PE")
            ax0.set_xlabel(r"Fuselage Location [m]")
            ax0.set_ylabel(f"x-location of the PE on the fuselage [m]")
            ax0.grid(True)
            ax0.set_xlim([0, l_fuselage])
            ax0.set_ylim([-1, l_fuselage * 1.05])
            ax0.legend()

            image2 = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\DUUC_diff_loc.png")

            ax1.imshow(image2, aspect='auto')
            ax1.set_axis_off()
            ax1.set_frame_on(False)

        if self.output_filepath:
            import os  # Make sure to import this at the top if not already
            filename = f"X_cg_plots_{plot_index}.png"
            full_path = os.path.join(self.output_filepath, filename)
            fig.savefig(full_path, dpi=300, bbox_inches='tight')
        canvas.draw()

    def plot_weight_distribution(self, w_vector1, w_vector2, prev_w_vector1, prev_w_vector2, plot_index, ax):
        labels = ['DUUC', 'ATR72-600']
        x = np.arange(len(labels))
        width = 0.35

        if plot_index == 0:  # Wing
            wing = np.array([w_vector1[1], w_vector2[1] + w_vector2[6] + w_vector2[7] + w_vector2[8]])
            ax.bar(x - width / 2, [w_vector1[1], 0], width, label='Wing', color='tab:blue')
            ax.bar(x - width / 2, [0, w_vector2[1]], width, bottom=[0, 0], color='tab:blue')
            ax.bar(x - width / 2, [0, w_vector2[6]], width, bottom=[0, w_vector2[1]], color='tab:red', label='Engine')
            ax.bar(x - width / 2, [0, w_vector2[7]], width, bottom=[0, w_vector2[1] + w_vector2[6]], color='tab:purple',
                   label='Fan')
            ax.bar(x - width / 2, [0, w_vector2[8]], width, bottom=[0, w_vector2[1] + w_vector2[6] + w_vector2[7]],
                   color='tab:orange', label='Nacelle')

            if prev_w_vector1 and prev_w_vector2:
                prev_wing = np.array(
                    [prev_w_vector1[1], prev_w_vector2[1] + prev_w_vector2[6] + prev_w_vector2[7] + prev_w_vector2[8]])
                ax.bar(x + width / 2, [prev_w_vector1[1], 0], width, color='tab:blue', alpha=0.5)
                ax.bar(x + width / 2, [0, prev_w_vector2[1]], width, bottom=[0, 0], color='tab:blue', alpha=0.5)
                ax.bar(x + width / 2, [0, prev_w_vector2[6]], width, bottom=[0, prev_w_vector2[1]], color='tab:red',
                       alpha=0.5)
                ax.bar(x + width / 2, [0, prev_w_vector2[7]], width, bottom=[0, prev_w_vector2[1] + prev_w_vector2[6]],
                       color='tab:purple', alpha=0.5)
                ax.bar(x + width / 2, [0, prev_w_vector2[8]], width,
                       bottom=[0, prev_w_vector2[1] + prev_w_vector2[6] + prev_w_vector2[7]], color='tab:orange',
                       alpha=0.5)

            ax.set_title('Wing')
            ax.set_ylim([0, 4000])
            ax.legend(loc="upper left")

        elif plot_index == 1:  # Empennage
            duuc_empennage = np.array([w_vector1[2:9]])
            atr_empennage = np.array([w_vector2[3:5]])

            # Plot DUUC Empennage components
            for i, component in enumerate(w_vector1[2:9]):
                ax.bar(x - width / 2, [component, 0], width, bottom=[sum(w_vector1[2:i + 2]), 0],
                       label=f'DUUC {["Duct", "Pylon", "Support", "Elevator", "Propeller Engine", "Nacelle", "Propeller Fan"][i]}',
                       color=plt.cm.tab20(i))

            # Plot ATR Empennage components
            for i, component in enumerate(w_vector2[3:5]):
                ax.bar(x - width / 2, [0, component], width, bottom=[0, sum(w_vector2[3:i + 3])],
                       label=f'ATR {["HT Tail", "VT Tail"][i]}', color=plt.cm.tab20(i + 7))

            if prev_w_vector1 and prev_w_vector2:
                prev_duuc_empennage = np.array([prev_w_vector1[2:9]])
                prev_atr_empennage = np.array([prev_w_vector2[3:5]])

                # Plot previous DUUC Empennage components
                for i, component in enumerate(prev_w_vector1[2:9]):
                    ax.bar(x + width / 2, [component, 0], width, bottom=[sum(prev_w_vector1[2:i + 2]), 0],
                           color=plt.cm.tab20(i), alpha=0.5)

                # Plot previous ATR Empennage components
                for i, component in enumerate(prev_w_vector2[3:5]):
                    ax.bar(x + width / 2, [0, component], width, bottom=[0, sum(prev_w_vector2[3:i + 3])],
                           color=plt.cm.tab20(i + 7), alpha=0.5)

            ax.set_title('Empennage')
            ax.set_ylim([0, 4000])
            ax.legend(loc="upper right")

        elif plot_index == 2:  # Aircraft Weight
            total_weight = np.array([sum(w_vector1), sum(w_vector2)])
            ax.bar(x - width / 2, total_weight, width, label='Total Weight', color='tab:blue')
            if prev_w_vector1 and prev_w_vector2:
                prev_total_weight = np.array([sum(prev_w_vector1), sum(prev_w_vector2)])
                ax.bar(x + width / 2, prev_total_weight, width, color='tab:blue', alpha=0.5)
            ax.set_title('Aircraft Weight')
            ax.set_ylim([0, 1.2 * max(total_weight)])
            ax.legend(loc="upper right")

        ax.set_xticks(x)
        ax.set_xticklabels(labels)

        ax.set_ylabel('Component mass [kg]')

    def update_cd0_empennage(self, plot_index, fig, ax, canvas):
        """Update CD0 empennage plots."""
        cd0_vector1 = self.calculation_results['CD0']['cd0_atr']
        cd0_vector2 = self.calculation_results['CD0']['cd0_duuc']

        # Reference values
        cd0_ref1 = np.array([1.702e-3, 2.065e-3])
        cd0_ref2 = np.array([49.2e-4, 2.8e-4, 0, 0])
        cd_totals = np.array([sum(cd0_vector1), sum(cd0_vector2)])
        cd_totals_ref = np.array([sum(cd0_ref1), sum(cd0_ref2)])

        label1 = ['Hor. Tail', 'Vert. Tail']
        label2 = ['Duct', 'Pylon', 'Support', 'Control']
        label3 = ['ATR72', 'DUUC']

        combined_labels = label1 + label2
        combined_cd0_vector = np.concatenate([cd0_vector1, cd0_vector2])
        combined_cd0_ref = np.concatenate([cd0_ref1, cd0_ref2])

        # Check if previous data exists
        has_previous_data = self.previous_calculation_results is not None

        # Adjust bar width and positions based on whether previous data exists
        if has_previous_data:
            width = 0.25  # Width for three sets of bars
        else:
            width = 0.35  # Width for two sets of bars - wider when only two sets

        if plot_index == 0:
            # Plot combined components of ATR72-600 and DUUC
            x = np.arange(len(combined_labels))

            if has_previous_data:
                # Three bars: current, previous, reference
                rects1 = ax.bar(x - width, combined_cd0_vector, width, label='Prediction model', color="tab:blue")

                prev_cd0_vector1 = self.previous_calculation_results['CD0']['cd0_atr']
                prev_cd0_vector2 = self.previous_calculation_results['CD0']['cd0_duuc']
                prev_combined_cd0_vector = np.concatenate([prev_cd0_vector1, prev_cd0_vector2])
                rects2 = ax.bar(x, prev_combined_cd0_vector, width, color="tab:blue", alpha=0.3,
                                label='Previous iteration')

                rects3 = ax.bar(x + width, combined_cd0_ref, width, label='Reference', color="tab:green")
            else:
                # Only two bars: current and reference, positioned symmetrically
                rects1 = ax.bar(x - width / 2, combined_cd0_vector, width, label='Prediction model', color="tab:blue")
                rects3 = ax.bar(x + width / 2, combined_cd0_ref, width, label='Reference', color="tab:green")

            ax.set_title('Comparison of Empennage Components')
            ax.set_ylabel('$C_{D0}$ [-]')
            ax.set_ylim([0, 0.01])
            ax.set_xticks(x)
            ax.set_xticklabels(combined_labels, rotation=45, ha='right')
            ax.legend(loc='upper left')

            # Set consistent padding around the plot
            fig.subplots_adjust(left=0.12, right=0.95, bottom=0.2, top=0.9)

        elif plot_index == 1:
            # Plot empennage totals
            x = np.arange(len(label3))

            if has_previous_data:
                # Three bars: current, previous, reference
                rects1 = ax.bar(x - width, cd_totals, width, label='Prediction model', color='tab:blue')

                prev_cd0_vector1 = self.previous_calculation_results['CD0']['cd0_atr']
                prev_cd0_vector2 = self.previous_calculation_results['CD0']['cd0_duuc']
                prev_cd_totals = np.array([sum(prev_cd0_vector1), sum(prev_cd0_vector2)])
                rects2 = ax.bar(x, prev_cd_totals, width, color='tab:blue', alpha=0.3,
                                label='Previous iteration')

                rects3 = ax.bar(x + width, cd_totals_ref, width, label='Reference', color="tab:green")
            else:
                # Only two bars: current and reference, positioned symmetrically
                rects1 = ax.bar(x - width / 2, cd_totals, width, label='Prediction model', color='tab:blue')
                rects3 = ax.bar(x + width / 2, cd_totals_ref, width, label='Reference', color="tab:green")

            ax.set_title('Empennage Totals')
            ax.set_ylabel('$C_{D0}$ [-]')
            ax.set_ylim([0, 0.01])
            ax.set_xticks(x)
            ax.set_xticklabels(label3, rotation=45, ha='right')
            ax.legend(loc='upper left')

            # Set consistent padding around the plot - matching the first plot
            fig.subplots_adjust(left=0.12, right=0.95, bottom=0.2, top=0.9)

        elif plot_index == 2:
            ax.clear()
            ax.set_facecolor("white")
            ax.set_axis_off()
            fig.patch.set_facecolor('white')
            reference_text = (r"$\bf{Reference\ Data\ Source:}$" "\n\n"
                              r"$\it{ATR72-600:}$" "Aircraft Design Studies \n "
                              "Based on the ATR 72, 2008\n\n"
                              r"$\it{DUUC:}$"
                              "Synthesis of an Aircraft Featuring a Ducted-Fan \n Propulsive Empennage, 2017\n")

            ax.text(0, 0.5, reference_text, horizontalalignment='left',
                    verticalalignment='center', transform=ax.transAxes, fontsize=12,
                    bbox=dict(facecolor='lightblue', edgecolor='black', alpha=0.7))

        elif plot_index == 2:
            ax.clear()
            ax.set_facecolor("white")
            ax.set_axis_off()
            fig.patch.set_facecolor('white')
            reference_text = (r"$\bf{Reference\ Data\ Source:}$" "\n\n"
                              r"$\it{ATR72-600:}$" "Aircraft Design Studies \n "
                              "Based on the ATR 72, 2008\n\n"
                              r"$\it{DUUC:}$"
                              "Synthesis of an Aircraft Featuring a Ducted-Fan \n Propulsive Empennage, 2017\n")

            ax.text(0, 0.5, reference_text, horizontalalignment='left',
                    verticalalignment='center', transform=ax.transAxes, fontsize=12,
                    bbox=dict(facecolor='lightblue', edgecolor='black', alpha=0.7))

        canvas.draw()
        if self.output_filepath:
            import os  # Make sure to import this at the top if not already
            filename = f"CD0_comparison_{plot_index}.png"
            full_path = os.path.join(self.output_filepath, filename)
            fig.savefig(full_path, dpi=300, bbox_inches='tight')

    def cleanup(self):
        """Clean up and disconnect all plots."""
        for component, plots in self.plots.items():
            for fig, ax, canvas in plots:
                canvas.setParent(None)
                plt.close(fig)
        self.plots = {}
        self.tab_widget.clear()

    def disconnect(self):
        """Disconnect any signal connections, if needed."""
        self.cleanup()
