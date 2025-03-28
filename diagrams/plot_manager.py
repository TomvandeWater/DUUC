import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt6.QtWidgets import QWidget, QGridLayout, QTabWidget
import numpy as np


class PlotManager(QWidget):
    def __init__(self, gui):
        super().__init__()
        self.gui = gui
        self.parameters = gui.parameters
        self.calculation_results = gui.calculation_results
        print(f"calculation results: {self.calculation_results}")
        self.plots = {}  # Dictionary to store plots for each component

        # Define all plot components
        self.components = [
            "Duct", "Pylon", "PE", "Support",
            "Weight", "X_cog", "Vtail", "Htail"
        ]

        # Initialize the layout
        self.layout = QGridLayout(self)

        # Create a tab widget
        self.tab_widget = QTabWidget()
        self.layout.addWidget(self.tab_widget, 0, 0)

        # Initialize tabs
        self.init_tabs()

    def init_tabs(self):
        """Initialize tabs for all components."""
        self.cleanup()  # Clean up old plots before initializing new ones

        for component in self.components:
            tab = QWidget()
            tab_layout = QGridLayout(tab)

            # Initialize three plots for each tab
            plots = []
            for j in range(3):
                fig, ax = plt.subplots(figsize=(4, 3))
                canvas = FigureCanvas(fig)
                tab_layout.addWidget(canvas, 0, j)
                plots.append((fig, ax, canvas))

            # Store the plots
            self.plots[component] = plots

            # Update plots based on the component
            self.update_plots(component)

            # Add the tab to the tab widget
            self.tab_widget.addTab(tab, component)

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
            elif component == "Weight":
                self.update_weight_plot(i, fig, ax, canvas)
            elif component in ["X_cog", "Vtail", "Htail"]:
                self.update_flight_mechanics_plot(component, i, fig, ax, canvas)

            ax.set_title(f"{component} - Plot {i + 1}")
            canvas.draw()

        print(f"Updated plots for {component}")  # Debug print

    def update_weight_plot(self, plot_index, fig, ax, canvas):
        """Update the weight plot."""
        w_vector_duuc = self.calculation_results['Weight']['w_vector_duuc']
        w_vector_atr = self.calculation_results['Weight']['w_vector_atr']
        self.plot_weight_distribution(w_vector_duuc, w_vector_atr, plot_index, ax)

    def update_flight_mechanics_plot(self, component, plot_index, fig, ax, canvas):
        """Update the plot for Flight Mechanics."""
        if component == "Weight":
            w_vector_duuc = self.calculation_results['Weight']['w_vector_duuc']
            w_vector_atr72_600 = self.calculation_results['Weight']['w_vector_atr']
            self.plot_weight_distribution(w_vector_duuc, w_vector_atr72_600, plot_index, ax)
        else:
            ax.text(0.5, 0.5, f"Flight Mechanics - {component} - Plot {plot_index}", ha='center', va='center')
        ax.set_title(f"{component} - Plot {plot_index}")
        canvas.draw()

    def plot_weight_distribution(self, w_vector1, w_vector2, plot_index, ax):
        labels = ['DUUC', 'ATR72-600']

        if plot_index == 0:  # Fuselage
            fuselage = np.array([w_vector1[0], w_vector2[0]])
            ax.bar(labels, fuselage, label='Fuselage', color='tab:blue')
            ax.set_title('Fuselage')
            ax.set_ylim([0, 6000])
        elif plot_index == 1:  # Wing
            wing = np.array([w_vector1[1], w_vector2[1]])
            ax.bar(labels, wing, label='Wing')
            ax.set_title('Wing')
            ax.set_ylim([0, 6000])
        elif plot_index == 2:  # Empennage
            htail = np.array([w_vector1[2], w_vector2[2]])
            vtail = np.array([w_vector1[3], w_vector2[3]])
            ax.bar(labels, htail, label='Horizontal tail')
            ax.bar(labels, vtail, label='Vertical tail', bottom=htail)
            ax.set_title('Empennage')
            ax.set_ylim([0, 2000])

        ax.legend(loc="lower right")
        ax.set_ylabel('Component mass [kg]')

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
