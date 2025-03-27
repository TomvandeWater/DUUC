import sys
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QTabWidget,
                             QGridLayout, QVBoxLayout, QHBoxLayout, QLabel,
                             QLineEdit, QGroupBox, QFormLayout, QComboBox, QMessageBox)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QIcon
from pyvistaqt import QtInteractor
from analysis_modules.aerodynamic import *
from visualize_aircraft import *
from plot_manager import PlotManager
from double_slide import DoubleSlider
import config
import data.atr_reference as ref
import ctypes


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Innovative Empennage Sizing')
        self.setWindowIcon(QIcon(r"C:\Users\tomva\pythonProject\DUUC\data\images\flame.ico"))

        # Store parameter values
        self.parameters = {
            'duct_diameter': config.duct_diameter, 'duct_chord': config.duct_chord, 'duct_profile': "NACA0016",
            'pylon_chord': config.pylon_chord, 'pylon_length': config.pylon_length,
            'support_chord': config.support_chord, 'support_length': config.support_length,
            'num_blades': config.n_blades, 'altitude': 7000,  'velocity': 10, 'alpha': 0.0, 'delta_e': 0, 'delta_r': 0,
            'power_condition': 'on', 'plot_group': 'Aerodynamics', 'propulsion_type': 'conventional',
            'cant_angle': config.cant_angle, 'pylon_profile': config.pylon_airfoil,
            'support_profile': config.support_airfoil, 'BEM1': 1,
            'RPM': config.rpm, 'propeller_diameter': 0.95 * config.duct_diameter, 'hub_diameter': config.hub_diameter,
            'propeller_airfoil': config.prop_airfoil, 'propeller_c_root': config.c_root, 'propeller_c_tip': config.c_tip,
            'wing_span': ref.b_w, 'wing_phi_qc': ref.phi_qc_w, 'wing_airfoil': ref.wing_airfoil, 'wing_tr': ref.tr_w,
            'wing_c_root': np.round(ref.c_root_w, 3),
            "fuselage_length": np.round((ref.l_tail+ref.l_cockpit+ref.l_cabin), 3),
            "fuselage_diameter": ref.diameter_fuselage, "fuselage_co_l": ref.l_cockpit, "fuselage_ca_l": ref.l_cabin,
            "fuselage_ta_l": ref.l_tail, "aircraft_n_pax": config.n_pax, "x_PE": np.round((12 + ref.l_tail), 3),
            "y_PE": 0, "z_PE": ref.diameter_fuselage, "nacelle_length": config.nacelle_length,
            "nacelle_diameter": config.nacelle_diameter, "hcv_span": config.control_vane_length, "hcv_chord": config.control_vane_chord,
            "vcv_span": config.control_vane_length, "vcv_chord": config.control_vane_chord, "cv_airfoil": config.control_vanes_airfoil,
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
        input_layout = QHBoxLayout(input_area)  # Changed to QHBoxLayout

        # Left side: Flight Conditions and Calculated Values
        left_side_layout = QVBoxLayout()

        # Flight Conditions Input (Top Right)
        flight_conditions_group = QGroupBox("Flight Conditions")
        flight_conditions_layout = QVBoxLayout()
        flight_conditions_group.setFixedWidth(700)
        self.create_flight_condition_inputs(flight_conditions_layout)

        self.plot_group_label = QLabel("Plot Group:")
        self.plot_group_selector = QComboBox()
        self.plot_group_selector.addItem("Aerodynamics")
        self.plot_group_selector.addItem("Flight Mechanics")
        self.plot_group_selector.setCurrentText(self.parameters['plot_group'])
        self.plot_group_selector.currentTextChanged.connect(self.update_plot_group)

        # Add to layout
        flight_conditions_layout.addWidget(self.plot_group_label)
        flight_conditions_layout.addWidget(self.plot_group_selector)

        flight_conditions_group.setLayout(flight_conditions_layout)
        left_side_layout.addWidget(flight_conditions_group)

        # Calculated Values Display
        self.calculated_values_group = QGroupBox("Calculated Values")
        calculated_values_layout = QFormLayout()
        self.calculated_values_group.setFixedWidth(700)
        self.mach_label = QLabel("Mach: ")
        self.density_label = QLabel("Density (kg/m^3): ")
        self.temperature_label = QLabel("Temperature (K): ")
        calculated_values_layout.addRow(self.mach_label)
        calculated_values_layout.addRow(self.density_label)
        calculated_values_layout.addRow(self.temperature_label)
        self.calculated_values_group.setLayout(calculated_values_layout)
        left_side_layout.addWidget(self.calculated_values_group)

        input_layout.addLayout(left_side_layout)  # Add left side layout to main input layout

        # Right Side: KPI and Calculation Status
        right_side_layout = QVBoxLayout()

        # Key Performance Indicators Display
        self.kpi_group = self.create_kpi_display()
        right_side_layout.addWidget(self.kpi_group)

        # Calculation Status Display
        self.calculation_status_group = self.create_calculation_status_display()
        right_side_layout.addWidget(self.calculation_status_group)

        input_layout.addLayout(right_side_layout)  # Add right side layout to main input layout

        # Parameter Input Tabs
        self.input_tabs = QTabWidget()
        left_side_layout.addWidget(self.input_tabs)  # Move this line up

        # Add parameter input tabs for components
        self.add_input_tab("Duct", self.create_duct_inputs)
        self.add_input_tab("Pylon", self.create_pylon_inputs)
        self.add_input_tab("Support", self.create_support_inputs)
        self.add_input_tab("Propeller", self.create_propeller_inputs)
        self.add_input_tab("Aircraft", self.create_aircraft_inputs)
        self.add_input_tab("Control Vane", self.create_control_inputs)
        self.add_input_tab("PE", self.create_pe_inputs)
        self.add_input_tab("BEM output", self.create_bem_output)

        self.layout.addWidget(input_area, 0, 1, 1, 1)  # Spans 1 rows, 1 column

        self.plot_manager = PlotManager(self.parameters, self.parameters['plot_group'])
        self.layout.addWidget(self.plot_manager, 1, 0, 1, 2)  # Span full width

        # Set row and column weights
        self.layout.setRowStretch(0, 1)  # Top row (3D view and input) gets more vertical space
        self.layout.setRowStretch(1, 1)  # Bottom row (plots)
        self.layout.setColumnStretch(0, 1)  # 3D view gets more horizontal space
        self.layout.setColumnStretch(1, 1)  # Input area

        # Timer for calculation status
        self.calculation_timer = QTimer(self)
        self.calculation_timer.timeout.connect(self.update_calculation_status_color)
        self.calculation_timer.start(100)  # Update every 100 ms

        self.is_calculating = False  # Flag to indicate calculation status
        self.is_plotting = False  # Flag to indicate plotting status

        # Initialize Status Color
        self.status_color = "green"
        self.status_color_plots = "green"
        self.update_status_color_labels()

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
            ref.b_w / 2,
            ref.c_root_w,
            11.25,
            11.5,
            25,
            15,
            self.parameters['x_PE'],
            self.parameters['y_PE'],
            self.parameters['z_PE'],
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
        plotter.enable_lightkit()
        visualize_aircraft(
            plotter,
            ref.diameter_fuselage,
            ref.b_w / 2,
            ref.c_root_w,
            11.25,
            11.5,
            25,
            15,
            0, 0, 0,
            "conventional",
            14,
        )
        plotter.reset_camera()

    def add_input_tab(self, tab_name, input_callback):
        tab_widget = QWidget()
        layout = QGridLayout(tab_widget)  # Changed layout to QGridLayout
        input_callback(layout)
        self.input_tabs.addTab(tab_widget, tab_name)

    def create_input_field(self, layout, row, col, label_text, param_name, data_type=float):
        label = QLabel(label_text)
        input_field = QLineEdit(str(self.parameters[param_name]))
        input_field.editingFinished.connect(lambda: self.update_parameter(param_name, input_field, data_type))
        layout.addWidget(label, row, col * 2)
        layout.addWidget(input_field, row, col * 2 + 1)

    def update_parameter(self, parameter, input_field, data_type=float):
        text = input_field.text()
        try:
            if data_type == float:
                value = float(text)
            elif data_type == int:
                value = int(text)
            else:
                value = text  # Treat as string

            self.parameters[parameter] = value
            self.start_calculation()  # Indicate calculation start
            QTimer.singleShot(500, self.finish_calculation)  # Simulate calculation time
            self.start_plotting()  # Indicate plotting start
            QTimer.singleShot(1000, self.finish_plotting)  # Simulate plotting time
            self.update_all_views()
            # self.update_all_plots()
            self.update_calculated_values()
            print(f"{parameter} updated to: {value}")  # For debugging
        except ValueError:
            QMessageBox.warning(self, "Input Error",
                                f"Invalid input for {parameter}. Please enter a valid {data_type.__name__}.")
            input_field.setText(str(self.parameters[parameter]))

    def create_flight_condition_inputs(self, layout):
        grid_layout = QGridLayout()

        # Existing inputs
        altitude_label = QLabel("Altitude [m]:")
        self.altitude_input = QLineEdit(str(self.parameters['altitude']))
        self.altitude_input.editingFinished.connect(
            lambda: self.update_flight_condition('altitude', self.altitude_input))

        velocity_label = QLabel("Velocity [m/s]:")
        self.velocity_input = QLineEdit(str(self.parameters['velocity']))
        self.velocity_input.editingFinished.connect(
            lambda: self.update_flight_condition('velocity', self.velocity_input))

        alpha_label = QLabel("Angle of Attack [deg]:")
        self.alpha_input = QLineEdit(str(self.parameters['alpha']))
        self.alpha_input.editingFinished.connect(
            lambda: self.update_flight_condition('alpha', self.alpha_input))

        # New dummy inputs
        ellev_deflect_label = QLabel("Elevator defl. [deg]:")
        self.ellev_deflect_input = QLineEdit(str(self.parameters['delta_e']))  # Dummy value
        self.ellev_deflect_input.editingFinished.connect(
            lambda: self.update_flight_condition('delta_e', self.ellev_deflect_input))

        rudder_deflect_label = QLabel("Rudder defl. [deg]:")
        self.rudder_deflect_input = QLineEdit(str(self.parameters['delta_r']))  # Dummy value
        self.rudder_deflect_input.editingFinished.connect(
            lambda: self.update_flight_condition('delta_r', self.rudder_deflect_input))

        # Add widgets to grid layout
        grid_layout.addWidget(altitude_label, 0, 0)
        grid_layout.addWidget(self.altitude_input, 0, 1)
        grid_layout.addWidget(velocity_label, 0, 2)
        grid_layout.addWidget(self.velocity_input, 0, 3)
        grid_layout.addWidget(alpha_label, 1, 0)
        grid_layout.addWidget(self.alpha_input, 1, 1)
        grid_layout.addWidget(ellev_deflect_label, 1, 2)
        grid_layout.addWidget(self.ellev_deflect_input, 1, 3)
        grid_layout.addWidget(rudder_deflect_label, 2, 0)
        grid_layout.addWidget(self.rudder_deflect_input, 2, 1)

        # Add power condition selector
        self.add_power_condition_selector(grid_layout)
        # Add propulsion type selector
        self.add_propulsion_type_selector(grid_layout)


        layout.addLayout(grid_layout)

        return self.altitude_input, self.velocity_input, self.alpha_input

    def add_power_condition_selector(self, layout):
        power_condition_label = QLabel("Power Condition:")
        self.power_condition_selector = QComboBox()
        self.power_condition_selector.addItems(["on", "off"])
        self.power_condition_selector.setCurrentText(self.parameters['power_condition'])
        self.power_condition_selector.currentTextChanged.connect(self.update_power_condition)

        layout.addWidget(power_condition_label, 2, 2)
        layout.addWidget(self.power_condition_selector, 2, 3)

    def update_power_condition(self, value):
        self.parameters['power_condition'] = value
        self.start_calculation()  # Indicate calculation start
        QTimer.singleShot(500, self.finish_calculation)  # Simulate calculation time
        self.start_plotting()  # Indicate plotting start
        QTimer.singleShot(1000, self.finish_plotting)  # Simulate plotting time
        # Add any additional logic needed when power condition changes
        print(f"Power condition updated to: {value}")

    def update_flight_condition(self, parameter, input_field):
        try:
            value = float(input_field.text())
            self.parameters[parameter] = value
            self.start_calculation()  # Indicate calculation start
            QTimer.singleShot(500, self.finish_calculation)  # Simulate calculation time
            self.start_plotting()  # Indicate plotting start
            QTimer.singleShot(1000, self.finish_plotting)  # Simulate plotting time
            self.update_calculated_values()  # Update calculated values when flight conditions change
            print(f"{parameter} updated to: {value}")  # For debugging
        except ValueError:
            # Handle the case where the input is not a valid float
            print(f"Invalid input for {parameter}")
            input_field.setText(str(self.parameters[parameter]))  # Revert to the original value

    def create_duct_inputs(self, layout):
        self.create_input_field(layout, 0, 0, "Duct Diameter:", "duct_diameter")
        self.create_input_field(layout, 0, 1, "Duct Chord:", "duct_chord")
        self.create_input_field(layout, 0, 2, "Profile:", "duct_profile", data_type=str)

    def create_pylon_inputs(self, layout):
        self.create_input_field(layout, 0, 0, "Pylon Chord:", "pylon_chord")
        self.create_input_field(layout, 0, 1, "Pylon Length:", "pylon_length")
        self.create_input_field(layout, 0, 2, "Cant angle:", "cant_angle")
        self.create_input_field(layout, 0, 3, "Airfoil:", "pylon_profile", data_type=str)

    def create_support_inputs(self, layout):
        self.create_input_field(layout, 0, 0, "Support Chord:", "support_chord")
        self.create_input_field(layout, 0, 1, "Support Length:", "support_length")
        self.create_input_field(layout, 0, 2, "Airfoil:", "support_profile", data_type=str)

    def create_propeller_inputs(self, layout):
        self.create_input_field(layout, 0, 0, "Propeller Diameter:", "propeller_diameter")
        self.create_input_field(layout, 0, 1, "Number of Blades:", "num_blades")
        self.create_input_field(layout, 0, 2, "RPM:", "RPM")
        self.create_input_field(layout, 0, 3, "Airfoil:", "propeller_airfoil", data_type=str)
        self.create_input_field(layout, 1, 0, "Hub Diameter:", "hub_diameter")
        self.create_input_field(layout, 1, 1, "Root chord:", "propeller_c_root")
        self.create_input_field(layout, 1, 2, "Tip Chord:", "propeller_c_tip")
        self.create_input_field(layout, 1, 3, "Nacelle diameter:", "nacelle_diameter")
        self.create_input_field(layout, 2, 0, "Nacelle length:", "nacelle_length")
        self.create_input_field(layout, 2, 1, "Propulsion type:", "propulsion_type")

    def create_control_inputs(self, layout):
        self.create_input_field(layout, 0, 0, "Hor. vane span:", "hcv_span")
        self.create_input_field(layout, 0, 1, "Hor. vane chord:", "hcv_chord")
        self.create_input_field(layout, 0, 2, "Airfoil:", "cv_airfoil", data_type=str)
        self.create_input_field(layout, 1, 0, "Ver. vane span:", "vcv_span")
        self.create_input_field(layout, 1, 1, "Ver. vane chord:", "vcv_chord")

    def create_aircraft_inputs(self, layout):
        self.create_input_field(layout, 0, 0, "Wing span:", "wing_span")
        self.create_input_field(layout, 0, 1, "Sweep:", "wing_phi_qc")
        self.create_input_field(layout, 0, 2, "NACA:", "wing_airfoil", data_type=str)
        self.create_input_field(layout, 0, 3, "Taper ratio:", "wing_tr")
        self.create_input_field(layout, 1, 0, "Root chord:", "wing_c_root")
        self.create_input_field(layout, 1, 1, "Fuselage length:", "fuselage_length")
        self.create_input_field(layout, 1, 2, "Fuselage diameter:", "fuselage_diameter")
        self.create_input_field(layout, 1, 3, "Cockpit length:", "fuselage_co_l")
        self.create_input_field(layout, 2, 0, "Cabin length:", "fuselage_ca_l")
        self.create_input_field(layout, 2, 1, "Tail length:", "fuselage_ta_l")
        self.create_input_field(layout, 2, 2, "Pax:", "aircraft_n_pax")

    def create_bem_output(self, layout):
        self.create_input_field(layout, 0, 0, "BEM1:", "BEM1")

    def create_pe_inputs(self, layout):
        self.create_input_field(layout, 0, 0, "X-position:", "x_PE")
        self.create_input_field(layout, 0, 1, "y-position:", "y_PE")
        self.create_input_field(layout, 0, 2, "z-position:", "z_PE")

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

    def create_kpi_display(self):
        kpi_group = QGroupBox("Key Performance Indicators")
        kpi_layout = QFormLayout()
        kpi_group.setFixedWidth(300)  # Set to desired width

        self.lift_label = QLabel("Lift [N]: ")
        self.drag_label = QLabel("Drag [N]: ")
        self.l_d_ratio_label = QLabel("L/D Ratio: ")
        self.thrust_label = QLabel("Thrust [N]: ")
        self.power_label = QLabel("Power [kW]: ")

        kpi_layout.addRow(self.lift_label)
        kpi_layout.addRow(self.drag_label)
        kpi_layout.addRow(self.l_d_ratio_label)
        kpi_layout.addRow(self.thrust_label)
        kpi_layout.addRow(self.power_label)

        kpi_group.setLayout(kpi_layout)
        return kpi_group

    def add_propulsion_type_selector(self, layout):
        propulsion_type_label = QLabel("Propulsion Type:")
        self.propulsion_type_selector = QComboBox()
        self.propulsion_type_selector.addItems(["conventional", "hybrid"])
        self.propulsion_type_selector.setCurrentText(self.parameters['propulsion_type'])
        self.propulsion_type_selector.currentTextChanged.connect(self.update_propulsion_type)

        layout.addWidget(propulsion_type_label, 3, 2)  # Adjust row and column as needed
        layout.addWidget(self.propulsion_type_selector, 3, 3)  # Adjust row and column as needed

    def update_propulsion_type(self, text):
        self.parameters['propulsion_type'] = text
        self.start_calculation()  # Indicate calculation start
        QTimer.singleShot(500, self.finish_calculation)  # Simulate calculation time
        self.start_plotting()  # Indicate plotting start
        QTimer.singleShot(1000, self.finish_plotting)  # Simulate plotting time
        self.update_all_views()
        # self.update_all_plots()
        self.update_calculated_values()
        print(f"Propulsion Type updated to: {text}")

    def create_calculation_status_display(self):
        status_group = QGroupBox("Calculation Status")
        status_group.setFixedHeight(75)  # Setting the height of the status box
        status_layout = QGridLayout()  # Using QGridLayout for two rows
        status_group.setFixedWidth(300)

        # Calculation Status
        self.status_label = QLabel("Status calculation:")
        self.status_color_label = QLabel()
        self.status_color_label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # Plots Status
        self.status_label_plots = QLabel("Status plots:")
        self.status_color_label_plots = QLabel()
        self.status_color_label_plots.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # Add to layout
        status_layout.addWidget(self.status_label, 0, 0)
        status_layout.addWidget(self.status_color_label, 0, 1)
        status_layout.addWidget(self.status_label_plots, 1, 0)
        status_layout.addWidget(self.status_color_label_plots, 1, 1)

        status_group.setLayout(status_layout)
        return status_group

    def update_status_color_labels(self):
        """Updates the color of the status labels."""
        style_sheet_calc = f"""
            QLabel {{
                background-color: {self.status_color};
                border-radius: 7px;  /* Make it round */
                min-width: 14px;      /* Set a minimum width */
                max-width: 14px;      /* Set a maximum width */
                min-height: 14px;     /* Set a minimum height */
                max-height: 14px;     /* Set a maximum height */
                margin-left: 5px;     /* Add some space to the left */
            }}
        """
        style_sheet_plots = f"""
            QLabel {{
                background-color: {self.status_color_plots};
                border-radius: 7px;  /* Make it round */
                min-width: 14px;      /* Set a minimum width */
                max-width: 14px;      /* Set a maximum width */
                min-height: 14px;     /* Set a minimum height */
                max-height: 14px;     /* Set a maximum height */
                margin-left: 5px;     /* Add some space to the left */
            }}
        """
        self.status_color_label.setStyleSheet(style_sheet_calc)
        self.status_color_label_plots.setStyleSheet(style_sheet_plots)

    def start_calculation(self):
        """Sets the status to orange and updates the GUI."""
        self.is_calculating = True
        self.status_color = "orange"
        self.update_status_color_labels()  # Update the color
        QApplication.processEvents()  # Force GUI update

    def finish_calculation(self):
        """Sets the status to green and updates the GUI."""
        self.is_calculating = False
        self.status_color = "green"
        self.update_status_color_labels()  # Update the color
        QApplication.processEvents()  # Force GUI update

    def start_plotting(self):
        """Sets the status to orange and updates the GUI."""
        self.is_plotting = True
        self.status_color_plots = "orange"
        self.update_status_color_labels()  # Update the color
        QApplication.processEvents()  # Force GUI update

    def finish_plotting(self):
        """Sets the status to green and updates the GUI."""
        self.is_plotting = False
        self.status_color_plots = "green"
        self.update_status_color_labels()  # Update the color
        QApplication.processEvents()  # Force GUI update

    def update_calculation_status_color(self):
        """Updates the color of the calculation status based on whether calculations are in progress."""
        if self.is_calculating:
            self.status_color = "orange"
        else:
            self.status_color = "green"

        if self.is_plotting:
            self.status_color_plots = "orange"
        else:
            self.status_color_plots = "green"
        self.update_status_color_labels()

    def update_plot_group(self, selected_plot_group):
        """Update the plot_group variable and re-initialize the PlotManager."""
        self.parameters['plot_group'] = selected_plot_group
        self.start_plotting()  # Indicate plotting start
        QTimer.singleShot(1000, self.finish_plotting)  # Simulate plotting time

        # 1. Remove the old plot_manager widget from the layout
        self.layout.removeWidget(self.plot_manager)

        # 2. Disconnect signals and prepare for deletion
        try:
            self.plot_manager.disconnect()  # Disconnect all signals (if possible)
        except AttributeError:
            pass  # If disconnect() method doesn't exist

        self.plot_manager.setParent(None)  # Disconnect from the parent (MainWindow)
        self.plot_manager.deleteLater()  # Mark the old plot_manager for deletion

        # 3. Create a new PlotManager instance
        self.plot_manager = PlotManager(self.parameters, selected_plot_group)

        # 4. Add the new plot_manager to the layout
        self.layout.addWidget(self.plot_manager, 1, 0, 1, 2)  # Span full width

    def update_calculated_values(self):
        altitude = self.parameters['altitude']
        velocity = self.parameters['velocity']
        alpha = self.parameters['alpha']

        temperature = air_density_isa(altitude)[1]
        density = air_density_isa(altitude)[0]
        mach = velocity / speed_of_sound(altitude)

        self.mach_label.setText(f"Mach [-]: {mach:>30.3f}")
        self.density_label.setText(f"Density [kg/m^3]: {density:>15.3f}")
        self.temperature_label.setText(f"Temperature [K]: {temperature:>20.3f}")

        # Dummy KPI calculations (not physically accurate)
        lift = density * velocity ** 2 * alpha * 10
        drag = density * velocity ** 2 * 0.05
        l_d_ratio = lift / drag if drag != 0 else 0
        thrust = drag * 1.2
        power = thrust * velocity / 1000  # Convert to kW

        self.lift_label.setText(f"Lift [N]: {lift:.2f}")
        self.drag_label.setText(f"Drag [N]: {drag:.2f}")
        self.l_d_ratio_label.setText(f"L/D Ratio: {l_d_ratio:.2f}")
        self.thrust_label.setText(f"Thrust [N]: {thrust:.2f}")
        self.power_label.setText(f"Power [kW]: {power:.2f}")

    def showEvent(self, event):
        super().showEvent(event)
        self.showMaximized()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    icon_path = r"C:\Users\tomva\pythonProject\DUUC\data\images\flame.ico"
    app.setWindowIcon(QIcon(icon_path))

    # Ensure Windows taskbar icon update
    if sys.platform == "win32":
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("DUUC_designer.myapp")

    window = MainWindow()
    window.show()
    sys.exit(app.exec())

