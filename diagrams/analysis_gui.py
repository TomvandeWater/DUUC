import sys
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QTabWidget,
                             QGridLayout, QVBoxLayout, QHBoxLayout, QLabel,
                             QLineEdit, QSizePolicy, QGroupBox, QFormLayout, QComboBox)
from PyQt6.QtCore import Qt
from pyvistaqt import QtInteractor
from analysis_modules.aerodynamic import *
from visualize_aircraft import *
from plot_manager import PlotManager
from double_slide import DoubleSlider


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Innovative Empennage Sizing')

        # Store parameter values
        self.parameters = {
            'duct_diameter': config.duct_diameter,
            'duct_chord': config.duct_chord,
            'pylon_chord': config.pylon_chord,
            'pylon_length': config.pylon_length,
            'support_chord': config.support_chord,
            'support_length': config.support_length,
            'propeller_diameter': 0.3 * config.duct_chord,
            'num_blades': config.n_blades,
            'altitude': 7000,  # Initial altitude
            'velocity': 10,  # Initial velocity
            'alpha': 0.0,
            'plot_group': 'Aerodynamics',
            'power_condition': 'On', # Default value for plot_group
        }

        # Main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.layout = QGridLayout(central_widget)  # Define layout as self.layout

        # 3D Visualization Area (Left)
        self.visualization_tabs = QTabWidget()
        self.layout.addWidget(self.visualization_tabs, 0, 0, 1, 1)  # Spans 1 rows, 1 column

        # Add tabs with different views
        self.add_visualization_tab("Master view", self.create_master_view)
        self.add_visualization_tab("Cross-section", self.create_cross_section)
        self.add_visualization_tab("Aircraft view", self.aircraft_view_DUUC)
        self.add_visualization_tab("Reference Aircraft", self.aircraft_view_ATR)

        # Input Area (Middle)
        input_area = QWidget()
        input_layout = QVBoxLayout(input_area)

        # Parameter Input Tabs
        self.input_tabs = QTabWidget()
        input_layout.addWidget(self.input_tabs)

        # Add parameter input tabs for components
        self.add_input_tab("Duct", self.create_duct_inputs)
        self.add_input_tab("Pylon", self.create_pylon_inputs)
        self.add_input_tab("Support", self.create_support_inputs)
        self.add_input_tab("Propeller", self.create_propeller_inputs)

        # Flight Conditions Input (Top Right)
        flight_conditions_group = QGroupBox("Flight Conditions")
        flight_conditions_layout = QFormLayout()
        flight_conditions_group.setFixedWidth(300)
        self.altitude_input, self.velocity_input, self.alpha_input = self.create_flight_condition_inputs(
            flight_conditions_layout)

        self.plot_group_label = QLabel("Plot Group:")
        self.plot_group_selector = QComboBox()
        self.plot_group_selector.addItem("Aerodynamics")
        self.plot_group_selector.addItem("Flight Mechanics")
        self.plot_group_selector.setCurrentText(self.parameters['plot_group'])
        self.plot_group_selector.currentTextChanged.connect(self.update_plot_group)

        # Add to layout
        flight_conditions_layout.addRow(self.plot_group_label, self.plot_group_selector)
        flight_conditions_group.setLayout(flight_conditions_layout)
        input_layout.addWidget(flight_conditions_group)

        # Calculated Values Display
        self.calculated_values_group = QGroupBox()
        calculated_values_layout = QFormLayout()
        self.calculated_values_group.setFixedWidth(300)  # Set to desired width
        self.mach_label = QLabel("Mach: ")
        self.density_label = QLabel("Density (kg/m^3): ")
        self.temperature_label = QLabel("Temperature (K): ")
        calculated_values_layout.addRow(self.mach_label)
        calculated_values_layout.addRow(self.density_label)
        calculated_values_layout.addRow(self.temperature_label)
        self.calculated_values_group.setLayout(calculated_values_layout)
        input_layout.addWidget(self.calculated_values_group)

        self.layout.addWidget(input_area, 0, 1, 1, 1)  # Spans 1 rows, 1 column

        self.plot_manager = PlotManager(self.parameters, self.parameters['plot_group'])
        self.layout.addWidget(self.plot_manager, 1, 0, 1, 2)  # Span full width

        # Set row and column weights
        self.layout.setRowStretch(0, 1)  # Top row (3D view and input) gets more vertical space
        self.layout.setRowStretch(1, 1)  # Bottom row (plots)
        self.layout.setColumnStretch(0, 1)  # 3D view gets more horizontal space
        self.layout.setColumnStretch(1, 1)  # Input area

        # self.update_all_plots()
        self.update_calculated_values()

    def add_visualization_tab(self, tab_name, visualization_callback):
        plotter = QtInteractor(self.visualization_tabs)
        self.visualization_tabs.addTab(plotter, tab_name)
        visualization_callback(plotter)

    def create_master_view(self, plotter):
        self.update_visualization(plotter)

    def create_cross_section(self, plotter):
        self.cross_section_view(plotter)
        plotter.camera_position = 'xy'
        plotter.camera.zoom(0)

    def update_visualization(self, plotter):
        plotter.clear()  # Clear existing geometry
        plotter.enable_lightkit()

        visualize_propulsive_empennage(
            plotter,
            self.parameters['pylon_chord'],
            self.parameters['pylon_length'],
            self.parameters['duct_chord'],
            self.parameters['duct_diameter'],
            self.parameters['support_chord'],
            self.parameters['support_length'],
            config.cant_angle,
            config.control_vane_chord,
            config.control_vane_length,
            0.25 * self.parameters['duct_chord'],
            0.95 * self.parameters['duct_chord'],
            0.40 * self.parameters['duct_chord'],
            0.30 * self.parameters['duct_chord']
        )
        plotter.reset_camera()

    def cross_section_view(self, plotter):
        plotter.clear()  # Clear existing geometry
        plotter.enable_lightkit()

        visualize_cross_section(
            plotter,
            self.parameters['pylon_chord'],
            self.parameters['pylon_length'],
            self.parameters['duct_chord'],
            self.parameters['duct_diameter'],
            self.parameters['support_chord'],
            self.parameters['support_length'],
            config.cant_angle,
            config.control_vane_chord,
            config.control_vane_length,
            0.25 * self.parameters['duct_chord'],
            0.95 * self.parameters['duct_chord'],
            0.40 * self.parameters['duct_chord'],
            0.30 * self.parameters['duct_chord']
        )
        plotter.reset_camera()

    def aircraft_view_DUUC(self, plotter):
        plotter.clear()
        plotter.enable_lightkit()

        visualize_aircraft(
            plotter,
            ref.diameter_fuselage,
            ref.b_w/2,
            ref.c_root_w,
            11.25,
            11.5,
            25,
            15,
            "DUUC",
            14,
            [self.parameters['pylon_chord'], self.parameters['pylon_length'], self.parameters['duct_chord'],
             self.parameters['duct_diameter'], self.parameters['support_chord'], self.parameters['support_length'],
             config.cant_angle, config.control_vane_chord, config.control_vane_length,
             0.25 * self.parameters['duct_chord'], 0.95 * self.parameters['duct_chord'],
             0.40 * self.parameters['duct_chord'], 0.30 * self.parameters['duct_chord']]
        )
        plotter.reset_camera()

    def aircraft_view_ATR(self, plotter):
        plotter.clear()
        visualize_aircraft(
            plotter,
            ref.diameter_fuselage,
            ref.b_w/2,
            ref.c_root_w,
            11.25,
            11.5,
            25,
            15,
            "conventional",
            14,
        )
        plotter.reset_camera()

    def add_input_tab(self, tab_name, input_callback):
        tab_widget = QWidget()
        layout = QVBoxLayout(tab_widget)
        input_callback(layout)
        self.input_tabs.addTab(tab_widget, tab_name)

    def create_slider_with_value(self, layout, label, param_name, min_val, max_val, decimals=2):
        label_widget = QLabel(label)
        slider = DoubleSlider(decimals=decimals, orientation=Qt.Orientation.Horizontal)
        slider.setMinimum(min_val)
        slider.setMaximum(max_val)
        slider.setValue(self.parameters[param_name])
        slider.setSingleStep(0.01)  # Set the step size for the slider

        value_display = QLineEdit(f"{self.parameters[param_name]:.{decimals}f}")
        value_display.setMaximumWidth(70)

        def update_value(value):
            self.parameters[param_name] = value
            value_display.setText(f"{value:.{decimals}f}")
            self.update_all_views()
            #self.update_all_plots()
            self.update_calculated_values()

        def update_slider():
            try:
                value = float(value_display.text())
                slider.setValue(value)
                self.parameters[param_name] = value
                self.update_all_views()
                #self.update_all_plots()
                self.update_calculated_values()
            except ValueError:
                pass

        slider.doubleValueChanged.connect(update_value)
        value_display.editingFinished.connect(update_slider)

        hlayout = QHBoxLayout()
        hlayout.addWidget(slider)
        hlayout.addWidget(value_display)

        form_layout = QFormLayout()
        form_layout.addRow(label_widget, hlayout)
        layout.addLayout(form_layout)

    def create_flight_condition_inputs(self, layout):
        altitude_label = QLabel("Altitude [m]:")
        altitude_input = QLineEdit(str(self.parameters['altitude']))
        altitude_input.editingFinished.connect(lambda: self.update_flight_condition('altitude', altitude_input))

        velocity_label = QLabel("Velocity [m/s]:")
        velocity_input = QLineEdit(str(self.parameters['velocity']))
        velocity_input.editingFinished.connect(lambda: self.update_flight_condition('velocity', velocity_input))

        alpha_label = QLabel("Angle of Attack [deg]:")
        alpha_input = QLineEdit(str(self.parameters['alpha']))
        alpha_input.editingFinished.connect(lambda: self.update_flight_condition('alpha', alpha_input))

        layout.addRow(altitude_label, altitude_input)
        layout.addRow(velocity_label, velocity_input)
        layout.addRow(alpha_label, alpha_input)

        return altitude_input, velocity_input, alpha_input

    def update_flight_condition(self, parameter, input_field):
        try:
            value = float(input_field.text())
            self.parameters[parameter] = value
            self.update_calculated_values()  # Update calculated values when flight conditions change
        except ValueError:
            # Handle the case where the input is not a valid float
            print(f"Invalid input for {parameter}")
            input_field.setText(str(self.parameters[parameter]))  # Revert to the original value

    def create_duct_inputs(self, layout):
        self.create_slider_with_value(layout, "Duct Diameter:", "duct_diameter", 0.0, 100.0)
        self.create_slider_with_value(layout, "Duct Chord:", "duct_chord", 0.0, 100.0)

    def create_pylon_inputs(self, layout):
        self.create_slider_with_value(layout, "Pylon Chord:", "pylon_chord", 0.0, 100.0)
        self.create_slider_with_value(layout, "Pylon Length:", "pylon_length", 0.0, 200.0)

    def create_support_inputs(self, layout):
        self.create_slider_with_value(layout, "Support Chord:", "support_chord", 0.0, 100.0)
        self.create_slider_with_value(layout, "Support Length:", "support_length", 0.0, 200.0)

    def create_propeller_inputs(self, layout):
        self.create_slider_with_value(layout, "Propeller Diameter:", "propeller_diameter", 0.0, 100.0)
        self.create_slider_with_value(layout, "Number of Blades:", "num_blades", 2, 8, decimals=0)

    def update_all_views(self):
        for i in range(self.visualization_tabs.count()):
            plotter = self.visualization_tabs.widget(i)
            if isinstance(plotter, QtInteractor):
                tab_name = self.visualization_tabs.tabText(i)

                if tab_name == "Master view":
                    self.create_master_view(plotter)
                elif tab_name == "Cross-section":
                    self.create_cross_section(plotter)
                elif tab_name == "Aircraft view":
                    self.aircraft_view_DUUC(plotter)
                elif tab_name == "Reference Aircraft":
                    self.aircraft_view_ATR(plotter)

    def update_plot_group(self, selected_plot_group):
        """Update the plot_group variable and re-initialize the PlotManager."""
        self.parameters['plot_group'] = selected_plot_group

        # 1. Remove the old plot_manager widget from the layout
        self.layout.removeWidget(self.plot_manager)

        # 2. Disconnect signals and prepare for deletion
        try:
            self.plot_manager.disconnect()  # Disconnect all signals (if possible)
        except AttributeError:
            pass  # If disconnect() method doesn't exist

        self.plot_manager.setParent(None)  # Disconnect from the parent (MainWindow)
        self.plot_manager.deleteLater()      # Mark the old plot_manager for deletion

        # 3. Create a new PlotManager instance
        self.plot_manager = PlotManager(self.parameters, selected_plot_group)

        # 4. Add the new plot_manager to the layout
        self.layout.addWidget(self.plot_manager, 1, 0, 1, 2)  # Span full width

    def update_calculated_values(self):
        altitude = self.parameters['altitude']
        velocity = self.parameters['velocity']

        temperature = air_density_isa(altitude)[1]
        density = air_density_isa(altitude)[0]
        mach = velocity / speed_of_sound(altitude)

        self.mach_label.setText(f"Mach [-]: {mach:>30.3f}")
        self.density_label.setText(f"Density [kg/m^3]: {density:>15.3f}")
        self.temperature_label.setText(f"Temperature [K]: {temperature:>20.3f}")

    def showEvent(self, event):
        super().showEvent(event)
        self.showMaximized()


if __name__ == "__main__":
    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()
    sys.exit(app.exec())
