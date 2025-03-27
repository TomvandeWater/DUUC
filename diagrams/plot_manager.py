import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from PyQt6.QtWidgets import QWidget, QGridLayout, QTabWidget


class PlotManager(QWidget):
    def __init__(self, parameters, plot_group):
        super().__init__()
        self.parameters = parameters
        self.plot_group = plot_group
        self.plots = {}  # Dictionary to store plots for each component

        # Define components for each plot group
        self.aerodynamics_components = ["Duct", "Pylon", "PE", "Support"]
        self.flight_mechanics_components = ["Weight", "X_cog", "Vtail", "Htail"]

        # Initialize the layout
        self.layout = QGridLayout(self)

        # Create a tab widget
        self.tab_widget = QTabWidget()
        self.layout.addWidget(self.tab_widget, 0, 0)

        # Initialize tabs based on the selected plot group
        self.init_tabs()

    def init_tabs(self):
        """Initialize tabs for the selected plot group."""
        self.cleanup()  # Clean up old plots before initializing new ones
        if self.plot_group == "Aerodynamics":
            components = self.aerodynamics_components
        elif self.plot_group == "Flight Mechanics":
            components = self.flight_mechanics_components

        for i, component in enumerate(components):
            tab = QWidget()
            tab_layout = QGridLayout(tab)

            # Initialize four plots for each tab
            plots = []
            for j in range(4):
                fig, ax = plt.subplots(figsize=(4, 3))
                canvas = FigureCanvas(fig)
                tab_layout.addWidget(canvas, 0, j)
                plots.append((fig, ax, canvas))

            # Store the plots
            self.plots[component] = plots

            # Update plots based on the component
            self.update_plots(component, plots)

            # Add the tab to the tab widget
            self.tab_widget.addTab(tab, component)

    def update_plots(self, component, plots):
        """Update the plots for the given component."""
        for i, (fig, ax, canvas) in enumerate(plots):
            ax.clear()
            if self.plot_group == "Aerodynamics":
                self.update_aerodynamics_plot(component, i, fig, ax, canvas)
            elif self.plot_group == "Flight Mechanics":
                self.update_flight_mechanics_plot(component, i, fig, ax, canvas)
            canvas.draw()  # Ensure the canvas is redrawn

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
            else:
                ax.text(0.5, 0.5, f"Plot {plot_index} for Duct", ha='center', va='center')  # Placeholder
        elif component == 'Pylon':
            if plot_index == 0:
                data = [self.parameters['pylon_chord'], self.parameters['pylon_length']]
                labels = ['Chord', 'Length']
                ax.bar(labels, data)
            else:
                ax.text(0.5, 0.5, f"Plot {plot_index} for Pylon", ha='center', va='center')  # Placeholder
        # Add other components as needed
        else:
            ax.text(0.5, 0.5, f"Plot {plot_index} for {component}", ha='center', va='center')  # Placeholder

        ax.set_title(f"{component} - Plot {plot_index}")
        canvas.draw()

    def update_flight_mechanics_plot(self, component, plot_index, fig, ax, canvas):
        """Update the plot for Flight Mechanics."""
        ax.text(0.5, 0.5, f"Flight Mechanics - {component} - Plot {plot_index}", ha='center', va='center')
        ax.set_title(f"{component} - Plot {plot_index}")
        canvas.draw()

    def cleanup(self):
        """Clean up and disconnect all plots."""
        for component, plots in self.plots.items():
            for fig, ax, canvas in plots:
                # Disconnect the canvas
                canvas.setParent(None)
                # Close the figure
                plt.close(fig)
        self.plots = {}  # Clear the stored plots
        self.tab_widget.clear()

    def disconnect(self):
        """Disconnect any signal connections, if needed."""
        # If there are specific signal connections in PlotManager, disconnect them here.
        # For example:
        # try:
        #     self.some_signal.disconnect(self.some_slot)
        # except TypeError:
        #     pass  # Signal was already disconnected or doesn't exist
        self.cleanup()
