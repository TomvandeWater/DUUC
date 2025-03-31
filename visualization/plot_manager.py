import matplotlib.ticker as ticker
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt6.QtWidgets import QWidget, QGridLayout, QTabWidget
from analysis_modules.factors import *
from PyQt6.QtGui import QIcon


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
            "Inflow",  "Geometry", "CD0", "Pylon", "PE", "Support",
            "Weight", "X_cog", "Vtail", "Htail", "requirements",
        ]

        # Initialize the layout
        self.layout = QGridLayout(self)

        # Create a tab widget
        self.tab_widget = QTabWidget()
        self.layout.addWidget(self.tab_widget, 0, 0)

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

            # Initialize three plots for each tab
            plots = []
            for j in range(3):
                if component == "X_cog" and j == 0:
                    fig, ax = plt.subplots(figsize=(12, 6))  # Wider figure for X_cog
                    canvas = FigureCanvas(fig)
                    tab_layout.addWidget(canvas, 0, 0, 1, 3)  # Span all three columns
                elif component == "X_cog" and j > 0:
                    # Skip creating additional plots for X_cog
                    continue
                else:
                    fig, ax = plt.subplots(figsize=(4, 3))  # Normal figure size for other plots
                    canvas = FigureCanvas(fig)
                    tab_layout.addWidget(canvas, 0, j)

                plots.append((fig, ax, canvas))

            # Store the plots
            self.plots[component] = plots

            # Update plots based on the component
            self.update_plots(component)

            # Add the tab to the tab widget
            self.tab_widget.addTab(tab, component)
            self.add_placeholder_icons(component)

    def update_aerodynamics_plot(self, component, plot_index, fig, ax, canvas):
        """Update the plot for Aerodynamics."""
        if component == 'Duct':
            if plot_index == 0:
                data = [self.parameters['duct_diameter'], self.parameters['duct_chord']]
                labels = ['Diameter', 'Chord']
                ax.bar(labels, data)
            elif plot_index == 1:
                data = [self.parameters['duct_diameter'], self.parameters['duct_chord']]
                labels = ['Diameter', 'Chord']
                ax.plot(labels, data)
        elif component == 'Pylon':
            if plot_index == 0:
                data = [self.parameters['pylon_chord'], self.parameters['pylon_length']]
                labels = ['Chord', 'Length']
                ax.bar(labels, data)
        ax.set_title(f"{component} - Plot {plot_index}")
        canvas.draw()

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
            if component in ["Duct", "Pylon", "PE", "Support"]:
                self.update_aerodynamics_plot(component, i, fig, ax, canvas)
            elif component == "Inflow":
                self.update_inflow_plot(i, fig, ax, canvas)
            elif component == "Weight":
                self.update_weight_plot(i, fig, ax, canvas)
            elif component == "CD0":
                self.update_cd0_empennage(i, fig, ax, canvas)
            elif component in ["X_cog", "Vtail", "Htail"]:
                self.update_flight_mechanics_plot(component, i, fig, ax, canvas)
            elif component in [ "Geometry"]:
                self.update_area_plot( i, fig, ax, canvas)

            canvas.draw()

        print(f"Updated plots for {component} tab")  # Debug print

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

        if component == "Vtail":
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

        canvas.draw()

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
        ax.set_ylim(0, 100)
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
        ax.set_ylim(0, 100)
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
        alpha = self.parameters['alpha']
        v_inf = self.parameters['velocity']

        v_prop = 0.94 * v_inf
        v_sup = 0.93 * v_inf
        v_cv = 0.97 * v_inf

        v_input2 = [v_inf, v_inf, v_prop, v_prop, v_sup, v_cv, (v_cv + v_inf) / 2, v_inf]
        a_input2 = [alpha, alpha/2, 0, 0, 0, 0, alpha/2, alpha]

        # Load background image
        img = plt.imread(r"C:\Users\tomva\pythonProject\DUUC\data\images\propulsive_empennage.png")

        if plot_index == 0:
            # Plot for inflow velocity
            ax.clear()
            ax.imshow(img, extent=[0.35, 6.00, 0.00, 200.00], aspect='auto')
            ax.plot(station, v_input, label=r"$V_{inflow}$ - power on", color="tab:blue", marker='o')
            ax.plot(station, v_input2, label=r"$V_{inflow}$ - power off", color="tab:blue", linestyle="dotted",
                    marker='o')
            ax.set_xlim(0, 6)
            ax.set_ylim(0, 200)
            ax.xaxis.set_major_locator(ticker.NullLocator())
            ax.set_ylabel(r'$V_{\infty}$ [m/s]')
            ax.set_title(r'Inflow Velocity per Component')
            ax.tick_params(axis='y')
            ax.grid(True)
            ax.legend(loc="upper left")

        elif plot_index == 1:
            # Plot for inflow angle
            ax.clear()
            ax.imshow(img, extent=[0.35, 6.00, -100, 100], aspect='auto')
            ax.plot(station, a_input, label=r"$\alpha$ - power on", color="tab:orange", marker='o')
            ax.plot(station, a_input2, label=r"$\alpha$ - power off", color="tab:orange", linestyle="dotted",
                    marker='o')
            ax.set_xlim(0, 6)
            ax.set_ylim(-100, 100)
            ax.xaxis.set_major_locator(ticker.NullLocator())
            ax.set_ylabel(r"$\alpha$ [deg]")
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
            ax.set_title(r'Area vs. chord (NACA0012)')

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
            ax.axvline(x=self.parameters["x_prop"], linestyle='--', color='black', alpha=0.5)
            ax.axhline(y=1.0, linestyle='--', color='tab:blue', alpha=0.5)

            # Add labels and title
            ax.set_xlabel(r'x/c [-]')
            ax.set_ylabel(r'$R_{i}$ / $R_{center}$ []')
            ax.set_title(r'Radius duct vs Radius available')

            # Add legend and grid
            ax.legend()
            ax.grid(True)

        elif plot_index == 2:
            # Empty plot for the third tab
            ax.clear()
            ax.set_axis_off()
            ax.set_facecolor('white')  # Set the background color to white
            fig.patch.set_facecolor('white')

        canvas.draw()

    def update_x_cog_plot(self, plot_index, fig, ax, canvas):
        """Update X_cog plot with horizontal position lines for both DUUC and ATR."""
        if plot_index == 0:
            # Get current and previous data

            current_duuc = self.calculation_results['X_cog']['x_cog_duuc']
            current_atr = self.calculation_results['X_cog']['x_cog_atr']
            prev_duuc = prev_atr = None
            if self.previous_calculation_results and 'X_cog' in self.previous_calculation_results:
                prev_duuc = self.previous_calculation_results['X_cog']['x_cog_duuc']
                prev_atr = self.previous_calculation_results['X_cog']['x_cog_atr']

            # Categories and styling
            categories = ['Fuselage CG', 'Wing CG', 'Overall CG', 'LEMAC']
            y_pos = np.arange(len(categories))
            offset = 0.2  # Vertical offset between DUUC and ATR lines

            ax.clear()
            ax.set_title('Center of Gravity Positions')
            # Plot current and previous values
            for i, (duuc_val, atr_val) in enumerate(zip(current_duuc, current_atr)):
                # DUUC current line
                ax.plot([0, duuc_val], [y_pos[i] + offset, y_pos[i] + offset],
                        color='tab:blue', linewidth=2, label='DUUC' if i == 0 else "")
                ax.plot(duuc_val, y_pos[i] + offset, marker='o', markersize=8, color='tab:blue')

                # ATR current line
                ax.plot([0, atr_val], [y_pos[i] - offset, y_pos[i] - offset],
                        color='tab:orange', linewidth=2, label='ATR' if i == 0 else "")
                ax.plot(atr_val, y_pos[i] - offset, marker='o', markersize=8, color='tab:orange')

                # Plot previous values if available
                if prev_duuc and prev_atr:
                    # DUUC previous line
                    ax.plot([0, prev_duuc[i]], [y_pos[i] + offset - 0.05, y_pos[i] + offset - 0.05],
                            color='tab:blue', alpha=0.3, linewidth=2)

                    # ATR previous line
                    ax.plot([0, prev_atr[i]], [y_pos[i] - offset + 0.05, y_pos[i] - offset + 0.05],
                            color='tab:orange', alpha=0.3, linewidth=2)

            # Configure plot appearance
            ax.set_yticks(y_pos)
            ax.set_yticklabels(categories)
            ax.set_xlabel('X Position [m]')

            ax.grid(True, axis='x')
            ax.legend(loc='lower right')

        elif plot_index in [1, 2]:
            # Empty plot for the second and third tabs
            ax.clear()
            ax.set_axis_off()  # Turn off the axis
            ax.set_facecolor('white')  # Set the background color to white
            fig.patch.set_facecolor('white')  # Set the figure background to white as well

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

        # Store previous data if available
        prev_cd0_vector1 = self.previous_calculation_results['CD0'][
            'cd0_atr'] if self.previous_calculation_results else None
        prev_cd0_vector2 = self.previous_calculation_results['CD0'][
            'cd0_duuc'] if self.previous_calculation_results else None

        cd0_ref1 = np.array([8.347e-4, 2.065e-3])
        cd0_ref2 = np.array([49.2e-4, 2.8e-4, 0, 0])
        cd_totals = np.array([sum(cd0_vector1), sum(cd0_vector2)])
        cd_totals_ref = np.array([sum(cd0_ref1), sum(cd0_ref2)])

        label1 = ['Hor. Tail', 'Vert. Tail']
        label2 = ['Duct', 'Pylon', 'Support', 'Control']
        label3 = ['ATR72-600', 'DUUC']

        combined_labels = label1 + label2
        combined_cd0_vector = np.concatenate([cd0_vector1, cd0_vector2])
        combined_cd0_ref = np.concatenate([cd0_ref1, cd0_ref2])

        width = 0.25  # Reduced width to accommodate previous iteration

        if plot_index == 0:
            # Plot combined components of ATR72-600 and DUUC
            x = np.arange(len(combined_labels))

            # Plot current data
            rects1 = ax.bar(x - width, combined_cd0_vector, width, label='Current', color="tab:blue")

            # Plot previous data with reduced opacity
            if prev_cd0_vector1 is not None and prev_cd0_vector2 is not None:
                prev_combined_cd0_vector = np.concatenate([prev_cd0_vector1, prev_cd0_vector2])
                rects2 = ax.bar(x, prev_combined_cd0_vector, width, color="tab:blue", alpha=0.3)

            # Plot reference data
            rects3 = ax.bar(x + width, combined_cd0_ref, width, label='Reference', color="tab:green")

            ax.set_title('Comparison of Empennage components')
            ax.set_ylabel('$C_{D0}$ [-]')
            ax.set_ylim([0, 0.01])
            ax.set_xticks(x)
            ax.set_xticklabels(combined_labels, rotation=45, ha='right')
            ax.legend(loc='upper left')

        elif plot_index == 1:
            # Plot empennage totals
            x = np.arange(len(label3))

            # Plot current data
            rects1 = ax.bar(x - width, cd_totals, width, label='Current', color='tab:blue')

            # Plot previous data with reduced opacity
            if prev_cd0_vector1 is not None and prev_cd0_vector2 is not None:
                prev_cd_totals = np.array([sum(prev_cd0_vector1), sum(prev_cd0_vector2)])
                rects2 = ax.bar(x, prev_cd_totals, width, color='tab:blue', alpha=0.3)

            # Plot reference data
            rects3 = ax.bar(x + width, cd_totals_ref, width, label='Reference', color="tab:green")

            ax.set_title('Empennage totals')
            ax.set_ylabel('$C_{D0}$ [-]')
            ax.set_xticks(x)
            ax.set_xticklabels(label3)
            ax.legend(loc='upper left')

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
