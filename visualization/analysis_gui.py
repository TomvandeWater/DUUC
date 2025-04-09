import time
from imports import *


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Innovative Empennage Sizing')
        self.setWindowIcon(QIcon(r"C:\Users\tomva\pythonProject\DUUC\data\images\flame.ico"))

        # Store parameter values
        self.parameters = {
            'duct_diameter': config.duct_diameter, 'duct_chord': config.duct_chord, 'duct_profile': config.duct_airfoil,

            'pylon_chord': config.pylon_chord, 'pylon_length': config.pylon_length, 'cant_angle': config.cant_angle,
            'pylon_profile': config.pylon_airfoil,

            'support_chord': config.support_chord, 'support_length': config.support_length,
            'support_profile': config.support_airfoil,

            'num_blades': config.n_blades, 'propeller_diameter': 0.95 * config.duct_diameter,
            'hub_diameter': config.hub_diameter, 'propeller_airfoil': config.prop_airfoil,
            'propeller_c_root': config.c_root, 'propeller_c_tip': config.c_tip,

            'BEM1': 41420, 'BEM2': 26482, 'BEM4': -1.44, 'BEM5': 0.889,
            'BEM6': 0.329, 'BEM7': 5, 'BEM8': 10, 'BEM3': 1820,

            'altitude': 7000,  'velocity': 128, 'alpha': 0.0, 'delta_e': 0, 'delta_r': 0,
            'power_condition': 'on', 'propulsion_type': 'conventional', 'RPM': config.rpm,
            "aircraft_n_pax": config.n_pax,

            'wing_span': ref.b_w, 'wing_phi_qc': ref.phi_qc_w, 'wing_airfoil': ref.wing_airfoil, 'wing_tr': ref.tr_w,
            'wing_c_root': np.round(ref.c_root_w, 3),

            "fuselage_length": np.round((ref.l_tail+ref.l_cockpit+ref.l_cabin), 3),
            "fuselage_diameter": ref.diameter_fuselage, "fuselage_co_l": ref.l_cockpit, "fuselage_ca_l": ref.l_cabin,
            "fuselage_ta_l": ref.l_tail,

            "nacelle_length": config.nacelle_length, "nacelle_diameter": config.nacelle_diameter,

            "hcv_span": config.control_vane_length, "hcv_chord": config.control_vane_chord,
            "vcv_span": config.control_vane_length, "vcv_chord": config.control_vane_chord,
            "cv_airfoil": config.control_vanes_airfoil,

            "l_v": 9.13, "x_prop": 0.3, "y_engine": ref.y_engine, "x_support": 0.5 * config.duct_chord,
            "x_control_vanes": 0.95 * config.duct_chord, "x_wing": 11.5, "x_pylon": 0.5 * config.duct_chord,
            "a_i_wing": 0, "a_i_duct": 0, "x_PE": 12, "y_PE": 0, "z_PE": ref.diameter_fuselage,
        }

        print("----- INITIALIZING XFOIL POLARS -----")
        self.update_xfoil_polars()  # initialized Xfoil polars for set parameters
        print("----- INITIALIZING POLARS COMPLETE -----")
        self.calculation_results = calculation_manager(self.parameters)
        self.previous_calculation_results = None

        # Main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Main vertical layout to hold taskbar + main content
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(1)

        # Add the taskbar at the top
        taskbar_widget = self.create_taskbar()
        main_layout.addWidget(taskbar_widget)

        # Create the main grid layout for content
        self.layout = QGridLayout()
        main_layout.addLayout(self.layout)

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

        flight_conditions_group.setLayout(flight_conditions_layout)
        left_side_layout.addWidget(flight_conditions_group)

        # --- Main Horizontal Layout to hold both sections ---
        calculated_and_terminal_layout = QHBoxLayout()

        # === Calculated Values Box ===
        self.calculated_values_group = QGroupBox("Calculated Values")
        self.calculated_values_group.setFixedWidth(350)  # Adjust width as needed

        # Main layout for both calculated values and BEM values
        main_layout = QVBoxLayout()

        # === Calculated Values Section ===
        values_layout = QGridLayout()

        self.mach_label = QLabel("Mach: ")
        self.density_label = QLabel("Density (kg/m^3): ")
        self.temperature_label = QLabel("Temperature (K): ")
        self.advance_label = QLabel("Advance ratio: ")

        # Add calculated values to grid layout
        values_layout.addWidget(self.mach_label, 0, 0)
        values_layout.addWidget(self.density_label, 0, 1)
        values_layout.addWidget(self.temperature_label, 1, 0)
        values_layout.addWidget(self.advance_label, 1, 1)

        # Add the calculated values grid layout to the main layout
        main_layout.addLayout(values_layout)

        # === BEM Values Section ===
        bem_header = QLabel("<b>BEM Module:</b>")  # Add header for BEM values
        main_layout.addWidget(bem_header)

        # Layout for BEM values
        bem_layout = QGridLayout()

        self.bem1_label = QLabel("BEM1: ")
        self.bem2_label = QLabel("BEM2: ")
        self.bem3_label = QLabel("BEM3: ")
        self.bem4_label = QLabel("BEM4: ")
        self.bem5_label = QLabel("BEM5: ")
        self.bem6_label = QLabel("BEM6: ")
        self.bem7_label = QLabel("BEM7: ")
        self.bem8_label = QLabel("BEM8: ")

        # Add BEM labels to the grid layout
        bem_layout.addWidget(self.bem1_label, 0, 0)
        bem_layout.addWidget(self.bem2_label, 0, 1)
        bem_layout.addWidget(self.bem3_label, 1, 0)
        bem_layout.addWidget(self.bem4_label, 1, 1)
        bem_layout.addWidget(self.bem5_label, 2, 0)
        bem_layout.addWidget(self.bem6_label, 2, 1)
        bem_layout.addWidget(self.bem7_label, 3, 0)
        bem_layout.addWidget(self.bem8_label, 3, 1)

        # Add the BEM layout to the main layout
        main_layout.addLayout(bem_layout)

        # Set the main layout for the group box (Calculated Values + BEM Values)
        self.calculated_values_group.setLayout(main_layout)

        # === Terminal Output Box ===
        self.terminal_output_group = QGroupBox("Terminal Prints")
        self.terminal_output_group.setFixedWidth(350)  # Adjust width as needed

        # Terminal layout
        terminal_layout = QVBoxLayout()

        self.command_display = QTextEdit()
        self.command_display.setReadOnly(True)
        self.command_display.setFixedHeight(150)
        terminal_layout.addWidget(self.command_display)

        self.terminal_output_group.setLayout(terminal_layout)

        # Add both sections (boxes) to the main layout
        calculated_and_terminal_layout.addWidget(self.calculated_values_group)
        calculated_and_terminal_layout.addWidget(self.terminal_output_group)

        # Add the combined layout to the left side layout
        left_side_layout.addLayout(calculated_and_terminal_layout)

        # Set up ConsoleOutput for terminal
        self.command_output = ConsoleOutput(self.command_display)
        sys.stdout = self.command_output

        # Add left side layout to the main layout
        input_layout.addLayout(left_side_layout)

        # Right Side: KPI and Calculation Status
        right_side_layout = QVBoxLayout()

        # Key Performance Indicators Display
        self.kpi_group = self.create_kpi_display()
        right_side_layout.addWidget(self.kpi_group)

        # Requirements Display
        self.requirements_group = self.create_requirements_display()
        right_side_layout.addWidget(self.requirements_group)

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

        self.layout.addWidget(input_area, 0, 1, 1, 1)  # Spans 1 rows, 1 column
        self.plot_manager = PlotManager(self)
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
        self.is_xfoiling = False  # Flag to indicate Xfoil status
        self.is_bemming = False  # Flag to indicate BEM status

        # Initialize Status Color
        self.status_color = "green"
        self.status_color_plots = "green"
        self.status_color_xfoil = "green"
        self.status_color_bem = "green"
        self.update_status_color_labels()

        # self.update_all_plots()
        self.update_calculated_values()

    def create_requirements_display(self):
        self.requirements_group = QGroupBox("Requirements")
        self.requirements_layout = QVBoxLayout()
        self.requirements_group.setLayout(self.requirements_layout)
        self.requirements_group.setFixedHeight(75)

        self.refresh_requirements_display()

        return self.requirements_group

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
            self.parameters['cant_angle'],
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
            self.parameters['cant_angle'],
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
            self.calculation_results["X_cog"]["x_cog_duuc"][3],
            self.calculation_results["X_cog"]["x_cog_duuc"][1],
            self.calculation_results["X_cog"]["x_cog_duuc"][0],
            self.calculation_results["X_cog"]["x_cog_duuc"][2],
            self.parameters['x_PE'] + self.calculation_results["X_cog"]["x_cog_duuc"][3],
            self.parameters['y_PE'],
            self.parameters['z_PE'],
            "DUUC",
            self.parameters['l_v'],
            0,
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
            self.calculation_results["X_cog"]["x_cog_atr"][3],
            self.calculation_results["X_cog"]["x_cog_atr"][1],
            self.calculation_results["X_cog"]["x_cog_atr"][0],
            self.calculation_results["X_cog"]["x_cog_atr"][2],
            0, 0, 0,
            "conventional",
            14, y_engine=self.parameters["y_engine"]
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

    def refresh_requirements_display(self):
        # Clear all existing widgets from the layout
        for i in reversed(range(self.requirements_layout.count())):
            item = self.requirements_layout.itemAt(i)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
                # If it's a layout, we need to clear and remove it
                self.clear_layout(item.layout())
                self.requirements_layout.removeItem(item)

        requirements = self.calculation_results.get('requirements', [])
        s_vert_requirement = requirements[0] if requirements else 0

        para = [
            (r"S-vert", np.round(s_vert_requirement, 2),
             2 * self.parameters.get("duct_diameter") * self.parameters.get("duct_chord")),
            (r"S-hor", 20, 2 * self.parameters.get("duct_diameter") * self.parameters.get("duct_chord")),
            ("Parameter 3", 30, 31),
            # Add more parameters as needed
        ]

        for name, requirement, actual in para:
            row_layout = QHBoxLayout()

            # Parameter name
            name_label = QLabel(name)
            row_layout.addWidget(name_label)

            # Requirement value
            req_label = QLabel(f"Req: {requirement}")
            row_layout.addWidget(req_label)

            # Actual value
            actual_label = QLabel(f"Actual: {actual}")
            row_layout.addWidget(actual_label)

            # Status indicator
            status_label = QLabel()
            status_label.setFixedSize(20, 20)
            status_label.setStyleSheet(
                f"""
                QLabel {{
                    background-color: {"green" if actual >= requirement else "red"};
                    border-radius: 7px;  /* Half of the width/height to make it circular */
                    min-width: 14px;
                    max-width: 14px;
                    min-height: 14px;
                    max-height: 14px;
                }}
                """
            )
            row_layout.addWidget(status_label)

            self.requirements_layout.addLayout(row_layout)

    def clear_layout(self, layout):
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.clear_layout(item.layout())

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
            xfoil_parameters = {'duct_chord', 'duct_profile', 'pylon_chord', 'support_chord',
                                'power_condition', 'pylon_profile', 'support_profile', 'BEM6', 'BEM7',
                                'hcv_chord', 'cv_airfoil'}
            if parameter in xfoil_parameters:
                self.start_xfoil()
                self.update_xfoil_polars()
                self.finish_xfoil()

            self.start_calculation()  # Indicate calculation start
            self.perform_calculation()
            self.finish_calculation()

            self.start_plotting()  # Indicate plotting start
            self.finish_plotting()
            self.update_all_views()

            self.update_calculated_values()
            self.refresh_requirements_display()
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

        bem_button = QPushButton("Run BEM analysis")
        bem_button.clicked.connect(self.run_bem_analysis)

        # Layout for buttons (side by side)
        button_layout = QHBoxLayout()
        button_layout.addWidget(bem_button)

        layout.addLayout(button_layout)  # Add buttons underneath everything

        return self.altitude_input, self.velocity_input, self.alpha_input

    def create_taskbar(self):
        taskbar = QWidget()
        taskbar.setFixedHeight(30)
        taskbar.setStyleSheet("background-color: #f0f0f0;")  # Light gray taskbar

        layout = QHBoxLayout()
        layout.setContentsMargins(10, 2, 10, 2)
        layout.setSpacing(5)

        save_button = QPushButton("💾 Save Configuration")
        save_button.clicked.connect(self.save_configuration_to_json)

        load_button = QPushButton("📂 Load Configuration")
        load_button.clicked.connect(self.load_configuration_from_json)

        layout.addWidget(save_button)
        layout.addWidget(load_button)
        layout.addStretch(1)  # Push buttons to the left

        taskbar.setLayout(layout)
        return taskbar

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
        self.perform_calculation()
        self.start_plotting()  # Indicate plotting start
        # Add any additional logic needed when power condition changes
        print(f"Power condition updated to: {value}")

    def update_flight_condition(self, parameter, input_field):
        try:
            value = float(input_field.text())
            self.parameters[parameter] = value
            xfoil_parameters = {'velocity', 'altitude'}
            if parameter in xfoil_parameters:
                self.start_xfoil()
                self.update_xfoil_polars()
                self.finish_xfoil()

            self.start_calculation()  # Indicate calculation start
            self.perform_calculation()
            self.finish_calculation()

            self.start_plotting()  # Indicate plotting start
            self.finish_plotting()

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
        self.create_input_field(layout, 0, 3, "X-location support:", 'x_support')

    def create_propeller_inputs(self, layout):
        self.create_input_field(layout, 0, 0, "Propeller Diameter:", "propeller_diameter")
        self.create_input_field(layout, 0, 1, "Number of Blades:", "num_blades")
        self.create_input_field(layout, 0, 2, "RPM:", "RPM")
        self.create_input_field(layout, 0, 3, "Airfoil:", "propeller_airfoil", data_type=str)
        self.create_input_field(layout, 1, 0, "Propeller loc.", "x_prop")
        self.create_input_field(layout, 1, 1, "Hub Diameter:", "hub_diameter")
        self.create_input_field(layout, 1, 2, "Root chord:", "propeller_c_root")
        self.create_input_field(layout, 1, 3, "Tip Chord:", "propeller_c_tip")
        self.create_input_field(layout, 2, 0, "Nacelle diameter:", "nacelle_diameter")
        self.create_input_field(layout, 2, 1, "Nacelle length:", "nacelle_length")
        self.create_input_field(layout, 2, 2, "Propulsion type:", "propulsion_type")

    def create_control_inputs(self, layout):
        self.create_input_field(layout, 0, 0, "Hor. vane span:", "hcv_span")
        self.create_input_field(layout, 0, 1, "Hor. vane chord:", "hcv_chord")
        self.create_input_field(layout, 0, 2, "Airfoil:", "cv_airfoil", data_type=str)
        self.create_input_field(layout, 1, 0, "Ver. vane span:", "vcv_span")
        self.create_input_field(layout, 1, 1, "Ver. vane chord:", "vcv_chord")
        self.create_input_field(layout, 1, 2, "X-location control:", 'x_control_vanes')

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
        self.create_input_field(layout, 2, 3, "Tail location:", "l_v")
        self.create_input_field(layout, 3, 0, "Y-engine ATR:", "y_engine")
        self.create_input_field(layout, 3, 1, "X-wing:", "x_wing")

    def create_bem_output(self, layout):
        self.create_input_field(layout, 0, 0, "BEM1:", "BEM1")
        self.create_input_field(layout, 0, 1, "BEM2:", "BEM2")
        self.create_input_field(layout, 0, 2, "BEM3:", "BEM3")
        self.create_input_field(layout, 0, 3, "BEM4:", "BEM4")
        self.create_input_field(layout, 1, 0, "BEM5:", "BEM5")
        self.create_input_field(layout, 1, 1, "BEM6:", "BEM6")
        self.create_input_field(layout, 1, 2, "BEM7:", "BEM7")
        self.create_input_field(layout, 1, 3, "BEM8:", "BEM8")
        self.create_input_field(layout, 2, 0, "BEM9:", "BEM9")

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
        self.update_calculated_values()
        print(f"Propulsion Type updated to: {text}")

    def create_calculation_status_display(self):
        status_group = QGroupBox("Prediction Model Status")
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

        # xfoil Status
        self.status_label_xfoil = QLabel("Status XFoil:")
        self.status_color_label_xfoil = QLabel()
        self.status_color_label_xfoil.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # BEM Status
        self.status_label_bem = QLabel("Status BEM:")
        self.status_color_label_bem = QLabel()
        self.status_color_label_bem.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # Add to layout
        status_layout.addWidget(self.status_label_xfoil, 0, 0)
        status_layout.addWidget(self.status_color_label_xfoil, 0, 1)
        status_layout.addWidget(self.status_label_bem, 0, 2)
        status_layout.addWidget(self.status_color_label_bem, 0, 3)
        status_layout.addWidget(self.status_label, 1, 0)
        status_layout.addWidget(self.status_color_label, 1, 1)
        status_layout.addWidget(self.status_label_plots, 1, 2)
        status_layout.addWidget(self.status_color_label_plots, 1, 3)

        status_group.setLayout(status_layout)
        return status_group

    def update_xfoil_polars(self):
        # -----                 NEW POLARS XFOIL                    ----- #
        altitude = self.parameters["altitude"]
        duct_chord = self.parameters['duct_chord']
        duct_profile = self.parameters['duct_profile']
        pylon_chord = self.parameters['pylon_chord']
        support_chord = self.parameters['support_chord']
        velocity = self.parameters['velocity']
        mach = velocity / speed_of_sound(altitude)
        power_condition = self.parameters['power_condition']
        pylon_profile = self.parameters['pylon_profile']
        support_profile = self.parameters['support_profile']
        BEM6 = self.parameters['BEM6']
        BEM7 = self.parameters['BEM7']
        hcv_chord = self.parameters['hcv_chord']
        cv_airfoil = self.parameters['cv_airfoil']

        density = air_density_isa(altitude)
        re_duct = reynolds(density, velocity, duct_chord)
        re_pylon = reynolds(density, velocity, pylon_chord)
        if power_condition == "off":
            v_sup = velocity
            v_control = velocity
        elif power_condition == "on":
            v_ax = 2 * abs(BEM6) + velocity
            v_tan = abs(BEM7)

            v_effective = np.sqrt(v_ax ** 2 + v_tan ** 2)
            v_sup = v_effective
            v_control = v_effective
        else:
            raise ValueError("Power Condition not properly specified")
        re_support = reynolds(density, v_sup, support_chord)
        re_control = reynolds(density, v_control, hcv_chord)
        return create_new_polars(profile_pylon=pylon_profile, profile_duct=duct_profile,
                                 profile_support=support_profile, profile_control=cv_airfoil, re_pylon=re_pylon,
                                 re_support=re_support, re_control=re_control, re_duct=re_duct, mach=mach)

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
        style_sheet_xfoil = f"""
                    QLabel {{
                        background-color: {self.status_color_xfoil};
                        border-radius: 7px;  /* Make it round */
                        min-width: 14px;      /* Set a minimum width */
                        max-width: 14px;      /* Set a maximum width */
                        min-height: 14px;     /* Set a minimum height */
                        max-height: 14px;     /* Set a maximum height */
                        margin-left: 5px;     /* Add some space to the left */
                    }}
                """
        style_sheet_bem = f"""
                    QLabel {{
                        background-color: {self.status_color_bem};
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
        self.status_color_label_xfoil.setStyleSheet(style_sheet_xfoil)
        self.status_color_label_bem.setStyleSheet(style_sheet_bem)

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

    def start_xfoil(self):
        """Sets the status to orange and updates the GUI."""
        self.is_xfoiling = True
        self.status_color_xfoil = "orange"
        self.update_status_color_labels()  # Update the color
        QApplication.processEvents()  # Force GUI update

    def finish_xfoil(self):
        """Sets the status to green and updates the GUI."""
        self.is_xfoiling = False
        self.status_color_xfoil = "green"
        self.update_status_color_labels()  # Update the color
        QApplication.processEvents()  # Force GUI update

    def start_bem(self):
        """Sets the status to orange and updates the GUI."""
        self.is_bemming = True
        self.status_color_bem = "orange"
        self.update_status_color_labels()  # Update the color
        QApplication.processEvents()  # Force GUI update

    def finish_bem(self):
        """Sets the status to green and updates the GUI."""
        self.is_bemming = False
        self.status_color_bem = "green"
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

        if self.is_xfoiling:
            self.status_color_xfoil = "orange"
        else:
            self.status_color_xfoil = "green"

        if self.is_bemming:
            self.status_color_bem = "orange"
        else:
            self.status_color_bem = "green"

        self.update_status_color_labels()

    def update_calculated_values(self):
        altitude = self.parameters['altitude']
        velocity = self.parameters['velocity']
        alpha = self.parameters['alpha']
        bem1 = self.parameters['BEM1']
        bem2 = self.parameters['BEM2']
        bem3 = self.parameters['BEM3']
        bem4 = self.parameters['BEM4']
        bem5 = self.parameters['BEM5']
        bem6 = self.parameters['BEM6']
        bem7 = self.parameters['BEM7']
        bem8 = self.parameters['BEM8']

        temperature = air_density_isa(altitude)[1]
        density = air_density_isa(altitude)[0]
        mach = velocity / speed_of_sound(altitude)
        advance_ratio = velocity / ((self.parameters["RPM"] / 60) * self.parameters["duct_diameter"])

        self.mach_label.setText(f"Mach [-]: {mach:>20.3f}")
        self.density_label.setText(f"Density [kg/m^3]: {density:>5.3f}")
        self.temperature_label.setText(f"Temperature [K]: {temperature:>8.1f}")
        self.advance_label.setText(f"Advance Ratio [-]: {advance_ratio:>6.1f}")
        self.bem1_label.setText(f"Bem1 [-]: {bem1:>6.1f}")
        self.bem2_label.setText(f"Bem2 [-]: {bem2:>6.1f}")
        self.bem3_label.setText(f"Bem3 [-]: {bem3:>6.1f}")
        self.bem4_label.setText(f"Bem4 [-]: {bem4:>6.1f}")
        self.bem5_label.setText(f"Bem5 [-]: {bem5:>6.1f}")
        self.bem6_label.setText(f"Bem6 [-]: {bem6:>6.1f}")
        self.bem7_label.setText(f"Bem7 [-]: {bem7:>6.1f}")
        self.bem8_label.setText(f"Bem8 [-]: {bem8:>6.1f}")

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

    def perform_calculation(self):
        try:
            # Store the previous calculation results
            if self.calculation_results is not None:
                self.previous_calculation_results = copy.deepcopy(self.calculation_results)

            self.calculation_results = calculation_manager(self.parameters)
            self.finish_calculation()

            self.start_plotting()
            # Update each component separately
            for component in self.calculation_results.keys():
                self.plot_manager.update_plots(component)

        except Exception as e:
            print(f"Calculation error: {e}")
            self.calculation_results = None
            self.previous_calculation_results = None

        finally:
            QTimer.singleShot(1000, self.finish_plotting)

    def save_configuration_to_json(self):
        # Prompt user for configuration name
        config_name, ok = QInputDialog.getText(self, "Save Configuration", "Enter configuration name:")

        if ok and config_name:
            # Define save path
            input_folder = r"C:\Users\tomva\pythonProject\DUUC\input"
            os.makedirs(input_folder, exist_ok=True)  # Ensure the folder exists

            json_path = os.path.join(input_folder, f"{config_name}.json")

            try:
                # Save as .json
                with open(json_path, "w") as json_file:
                    json.dump(self.parameters, json_file, indent=4)

                QMessageBox.information(self, "Success", f"Configuration saved to:\n{json_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save configuration:\n{e}")

    def load_configuration_from_json(self):
        input_folder = r"C:\Users\tomva\pythonProject\DUUC\input"
        dialog_result = QFileDialog.getOpenFileName(self, "Load Configuration", input_folder, "JSON Files (*.json)")

        if dialog_result and dialog_result[0]:
            file_path = dialog_result[0]
            try:
                # Read parameters exactly as-is from JSON
                with open(file_path, 'r') as file:
                    file_params = json.load(file)

                # Update self.parameters with raw values
                self.parameters.update(file_params)

                # Refresh GUI and calculations
                self.refresh_all_input_tabs()
                self.start_calculation()
                self.perform_calculation()
                self.finish_calculation()
                self.update_all_views()
                self.update_calculated_values()
                self.refresh_requirements_display()

                QMessageBox.information(self, "Success", f"Configuration loaded from:\n{file_path}")

            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load configuration:\n{str(e)}")

    def run_bem_analysis(self):
        self.start_bem()
        try:
            print(colored("Matlab is starting.....", "blue"))
            import matlab.engine
            bem_matlab_engine = matlab.engine.start_matlab()
            bem_matlab_engine.cd(r'C:\Users\tomva\pythonProject\DUUC\analysis_modules\BEM')
            print("BEM model has started")
            altitude = self.parameters['altitude']
            velocity = self.parameters['velocity']
            advance_ratio = velocity / ((self.parameters["RPM"] / 60) * self.parameters["duct_diameter"])

            temperature = air_density_isa(altitude)[1]
            density = air_density_isa(altitude)[0]

            t_out, q_out, n_out, tc, cp, ct, va, vt = bem_matlab_engine.BEM2(self.parameters["num_blades"],
                                                                             self.parameters["propeller_diameter"],
                                                                             25, 0, velocity, 0, advance_ratio,
                                                                             density, temperature,
                                                                             self.parameters["propeller_airfoil"],
                                                                             nargout=8)
            # linking the output variables back to the gui outputs
            self.parameters["BEM1"] = t_out
            self.parameters["BEM2"] = q_out
            self.parameters["BEM3"] = n_out
            self.parameters["BEM4"] = tc
            self.parameters["BEM5"] = cp
            self.parameters["BEM6"] = ct
            self.parameters["BEM7"] = va
            self.parameters["BEM8"] = vt

            print(colored("Success: BEM model completed", "green"))
            bem_matlab_engine.quit()
            print(colored("Matlab closes.....", "blue"))
        except Exception as e:
            import traceback
            error_message = traceback.format_exc()
            print("An error occurred during BEM analysis:\n", error_message)
            QMessageBox.critical(self, "BEM Analysis Error", f"An error occurred:\n{str(e)}")

        self.finish_bem()
        self.start_xfoil()
        self.update_xfoil_polars()
        self.finish_xfoil()

        self.start_calculation()  # Indicate calculation start
        self.perform_calculation()
        self.finish_calculation()

        self.start_plotting()  # Indicate plotting start
        self.finish_plotting()
        self.update_all_views()

        self.update_calculated_values()
        self.refresh_requirements_display()

    def refresh_all_input_tabs(self):
        """Refreshes all input tabs to reflect the updated parameters."""
        for i in range(self.input_tabs.count()):
            tab_widget = self.input_tabs.widget(i)
            if tab_widget:
                layout = tab_widget.layout()
                if layout:
                    # Clear existing widgets from the layout
                    for i in reversed(range(layout.count())):
                        item = layout.itemAt(i)
                        widget = item.widget()
                        if widget is not None:
                            widget.deleteLater()
                    # Re-create the input fields in the tab
                    tab_name = self.input_tabs.tabText(i)
                    if tab_name == "Duct":
                        self.create_duct_inputs(layout)
                    elif tab_name == "Pylon":
                        self.create_pylon_inputs(layout)
                    elif tab_name == "Support":
                        self.create_support_inputs(layout)
                    elif tab_name == "Propeller":
                        self.create_propeller_inputs(layout)
                    elif tab_name == "Aircraft":
                        self.create_aircraft_inputs(layout)
                    elif tab_name == "Control Vane":
                        self.create_control_inputs(layout)
                    elif tab_name == "PE":
                        self.create_pe_inputs(layout)
                    elif tab_name == "BEM output":
                        self.create_bem_output(layout)

    def showEvent(self, event):
        super().showEvent(event)
        self.showMaximized()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    icon_path = r"C:\Users\tomva\pythonProject\DUUC\data\images\flame.ico"

    # Set taskbar icon via AUMID
    if sys.platform == "win32":
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("DUUC_designer2.myapp")

    # Set window icon
    app.setWindowIcon(QIcon(icon_path))

    # Create main window
    window = MainWindow()
    window.move(screen_param()[0], screen_param()[1])
    window.setWindowIcon(QIcon(icon_path))  # Set window icon explicitly
    window.show()

    # Start application event loop
    sys.exit(app.exec())

