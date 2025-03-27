import sys
from qtpy import QtWidgets
import pyvista as pv
from pyvistaqt import QtInteractor
from visualize_propulsive_empennage import *


# Main application class
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, width=1600, height=900):
        super().__init__()
        self.setWindowTitle("Propulsive Empennage")
        self.resize(width, height)  # Set the window size

        # Create a tab widget
        self.tabs = QtWidgets.QTabWidget()
        self.setCentralWidget(self.tabs)

        # Add tabs with PyVista plots
        self.add_tab("Master view", lambda plotter: visualize_propulsive_empennage(
            plotter,
            config.pylon_chord,
            config.pylon_length,
            config.duct_chord,
            config.duct_diameter,
            config.support_chord,
            config.support_length,
            config.cant_angle,
            config.control_vane_chord,
            config.control_vane_length,
            0.25 * config.duct_chord,
            0.95 * config.duct_chord,
            0.40 * config.duct_chord,
            0.30 * config.duct_chord
        ))

        self.add_tab("Cross-section", lambda plotter: visualize_cross_section(
            plotter,
            config.pylon_chord,
            config.pylon_length,
            config.duct_chord,
            config.duct_diameter,
            config.support_chord,
            config.support_length,
            config.cant_angle,
            config.control_vane_chord,
            config.control_vane_length,
            0.25 * config.duct_chord,
            0.95 * config.duct_chord,
            0.40 * config.duct_chord,
            0.30 * config.duct_chord
        ))

        self.add_tab("Parameters", lambda plotter: get_parameters(
            plotter,
            config.pylon_chord,
            config.pylon_length,
            config.duct_chord,
            config.duct_diameter,
            config.support_chord,
            config.support_length,
            config.cant_angle,
            config.control_vane_chord,
            config.control_vane_length,
            0.25 * config.duct_chord,
            0.95 * config.duct_chord,
            0.40 * config.duct_chord,
            0.30 * config.duct_chord
        ))

    def add_tab(self, title, plot_function):
        # Create a QWidget for the tab
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)

        # Add a PyVista QtInteractor to the tab
        plotter = QtInteractor(tab)
        layout.addWidget(plotter.interactor)

        # Call the provided plotting function
        plot_function(plotter)

        # Add the tab to the QTabWidget
        self.tabs.addTab(tab, title)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


