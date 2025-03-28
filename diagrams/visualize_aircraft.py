import pyvista as pv
import numpy as np
import data.atr_reference as ref
from visualize_propulsive_empennage import *
import config


def visualize_propulsive_empennage_ac(c_pylon, b_pylon, c_duct, d_duct, c_support, b_support, cant, c_cv, b_cv, LE_pylon,
                                      LE_control, LE_support, x_prop, origin=(0, 0, 0)):

    primary = cmyk_to_rgb(1, 0, 0, 0)
    off_white = cmyk_to_rgb(0, 0, 0.02, 0.02)
    file_path_0012 = os.path.join(r"C:\Users\tomva\pythonProject\DUUC\data\airfoil_coordinates", "Naca0012" + ".txt")
    cant_rad = np.radians(cant)

    # Displace the origin by adding the origin offsets
    origin_x, origin_y, origin_z = origin

    # Pylon geometry
    pylon = plot_straight_wing(file_path_0012, span=np.cos(cant_rad) * b_pylon, chord=c_pylon, dihedral_deg=cant,
                               start_point=(LE_pylon + origin_x, origin_y, origin_z))

    x_pylon_end = LE_support + origin_x
    y_pylon_end = origin_y + np.cos(cant_rad) * b_pylon
    z_pylon_end = origin_z + np.sin(cant_rad) * b_pylon  # Removed extra origin_z addition

    # Support geometry
    support = plot_straight_wing(file_path_0012, span=np.cos(cant_rad) * b_support, chord=c_support, dihedral_deg=cant,
                                 start_point=(x_pylon_end, y_pylon_end, z_pylon_end))

    # Duct placement
    x_duct = origin_x
    y_duct = y_pylon_end + np.cos(cant_rad) * (0.5 * d_duct)
    z_duct = z_pylon_end + np.sin(cant_rad) * (0.5 * d_duct)  # Removed redundant origin_z

    duct = visualize_annular_wing(file_path_0012, chord=c_duct, diameter=d_duct,
                                  start_point=(x_duct, y_duct, z_duct), num_sections=100)

    # Apply logo texture to the duct
    logo_width = c_duct * 0.75
    logo_height = c_duct * 0.4
    logo_center = (x_duct + 0.15 * c_duct, y_duct + d_duct / 2, z_duct - 0.5 * logo_height)

    duct_with_logo, logo_texture = apply_single_logo_texture(
        duct,
        r"C:\Users\tomva\pythonProject\DUUC\data\images\TU_Delft.png",
        center=logo_center,
        size_u=logo_width,
        size_v=logo_height,
        direction_u=(1, 0, 0),
        direction_v=(0, 0, 1)
    )

    # Horizontal control vane
    x_control_h = LE_control + origin_x
    y_control_h = y_pylon_end + np.cos(cant_rad) * b_support / 2 - 0.5 * d_duct
    z_control_h = z_pylon_end + np.sin(cant_rad) * b_support / 2  # Removed redundant origin_z

    h_control = plot_straight_wing(file_path_0012, span=b_cv * 2, chord=c_cv, dihedral_deg=0,
                                   start_point=(x_control_h, y_control_h, z_control_h))

    # Vertical control vane
    x_control_v = LE_control + origin_x
    y_control_v = y_pylon_end + np.cos(cant_rad) * b_support / 2
    z_control_v = z_pylon_end + np.sin(cant_rad) * b_support / 2 - 0.5 * d_duct  # Removed redundant origin_z

    v_control = plot_vertical_wing(file_path_0012, span=b_cv * 2, chord=c_cv, sweep_deg=0.0,
                                   start_point=(x_control_v, y_control_v, z_control_v))

    # Engine components
    x_engine = x_prop + origin_x
    spinner = plot_half_sphere(config.hub_diameter, (x_engine, y_duct, z_duct))

    propeller = create_propeller(file_path_0012, span=d_duct / 2, center_point=(x_engine - config.c_root, y_duct, z_duct),
                                 num_wings=config.n_blades, chord=config.c_root)

    nacelle = create_engine_nacelle(diameter=config.hub_diameter, length=config.nacelle_length,
                                    cone_length=0.5, center_point=(x_engine, y_duct, z_duct), resolution=100)

    combined = pylon + support + h_control + v_control + duct_with_logo + spinner + propeller + nacelle
    return combined


def plot_vertical_wing_t(file_path, span=1.0, root_chord=1.0, tip_chord=0.5, sweep_deg=0.0, start_point=(0, 0, 0),
                         normalized=True):
    airfoil_coords = np.loadtxt(file_path)
    x_raw, y_raw = airfoil_coords[:, 0], airfoil_coords[:, 1]

    # Root profile
    x_root = x_raw * root_chord
    y_root = y_raw * root_chord if normalized else y_raw
    root_profile = np.column_stack((x_root, y_root, np.zeros_like(x_root)))

    # Tip profile with taper
    x_tip = x_raw * tip_chord
    y_tip = y_raw * tip_chord if normalized else y_raw

    # Sweep calculation (applied in X direction)
    sweep_rad = np.radians(sweep_deg)
    tip_x_offset = span * np.tan(sweep_rad)  # Sweep offset applied in X direction
    tip_z_offset = span  # The tip is at the end of the span in Z direction

    tip_profile = np.column_stack((x_tip + tip_x_offset, y_tip, np.full_like(x_tip, tip_z_offset)))

    # Translate to start point
    root_profile += np.array(start_point)
    tip_profile += np.array(start_point)

    # Combine points
    all_points = np.vstack((root_profile, tip_profile))
    n_points = len(root_profile)

    # Quad faces
    faces = []
    for i in range(n_points - 1):
        p1, p2, p3, p4 = i, i + n_points, i + 1 + n_points, i + 1
        faces.append([4, p1, p2, p3, p4])

    # Root cap
    root_center = np.mean(root_profile, axis=0)
    root_center_idx = len(all_points)
    all_points = np.vstack((all_points, root_center))
    for i in range(n_points - 1):
        faces.append([3, root_center_idx, i, i + 1])

    # Tip cap
    tip_center = np.mean(tip_profile, axis=0)
    tip_center_idx = len(all_points)
    all_points = np.vstack((all_points, tip_center))
    for i in range(n_points - 1):
        faces.append([3, tip_center_idx, i + 1 + n_points, i + n_points])

    faces = np.hstack(faces)
    wing_mesh = pv.PolyData(all_points, faces)

    return wing_mesh


def generate_cylinder(centers, radii, num_points=50):
    # Generate cross-section points at different x-coordinates in the y-z plane
    points = []
    for (x, y_center, z_center), r in zip(centers, radii):
        theta = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
        y = y_center + r * np.cos(theta)
        z = z_center + r * np.sin(theta)
        x_arr = np.full_like(y, x)
        points.extend(np.column_stack([x_arr, y, z]))

    # Convert to numpy array
    points = np.array(points)

    # Create faces (connecting circles into a closed surface)
    cells = []
    for i in range(len(centers) - 1):
        for j in range(num_points):
            next_j = (j + 1) % num_points
            lower = i * num_points + j
            upper = (i + 1) * num_points + j
            cells.append([4, lower, upper, upper + 1 if next_j != 0 else (i + 1) * num_points, lower + 1 if next_j != 0 else i * num_points])

    # Create PolyData mesh
    faces = np.hstack(cells).astype(np.int32)
    cylinder = pv.PolyData(points, faces)

    return cylinder


def visualize_aircraft(plotter, d_fus, semi_b_wing, c_wing, x_lemac, x_cog_w, x_cog_f, x_cog, x_duuc, y_duuc, z_duuc, aircraft_type: str, lv,
                       pe_input=None):
    off_white = cmyk_to_rgb(0, 0, 0.02, 0.02)
    r_fus = d_fus / 2
    c_cockpit = [(0, 0, 0.4*r_fus), (0.3*ref.l_cockpit, 0, 0.5*r_fus), (0.5*ref.l_cockpit, 0, 0.7*r_fus), (ref.l_cockpit, 0, r_fus)]

    radi_cockpit = [0.4*r_fus, 0.5 * r_fus, 0.7 * r_fus, r_fus]
    cockpit = generate_cylinder(c_cockpit, radi_cockpit)
    c_fuselage = [(ref.l_cockpit, 0, r_fus), (ref.l_cockpit + ref.l_cab, 0, r_fus)]
    radi_fus = [r_fus, r_fus]

    fuselage = generate_cylinder(c_fuselage, radi_fus)

    c_tail = [(ref.l_cockpit + ref.l_cab, 0, r_fus), (ref.l_cockpit + ref.l_cab + ref.l_tail, 0, 1.7 * r_fus)]
    radi_tail = [r_fus, 0.3 * r_fus]

    tail = generate_cylinder(c_tail, radi_tail)
    raidome = plot_half_sphere(0.8*r_fus, (0, 0, 0.4*r_fus))
    body = tail + raidome + cockpit + fuselage

    file_path_0012 = os.path.join(r"C:\Users\tomva\pythonProject\DUUC\data\airfoil_coordinates", "Naca0012" + ".txt")
    wing1 = plot_straight_wing(file_path_0012, span=semi_b_wing, chord=c_wing, dihedral_deg=0,
                               start_point=(x_lemac, 0, 1.90*r_fus))
    wing2 = plot_straight_wing(file_path_0012, span=-semi_b_wing, chord=c_wing, dihedral_deg=0,
                               start_point=(x_lemac, 0, 1.90*r_fus))
    wing_comb = wing1 + wing2

    # ---- create Cog Lines -----#
    z_min, z_max = 0, 10  # Adjust as needed
    line = pv.Line((x_cog_w, 0, z_min), (x_cog_w, 0, z_max))
    line2 = pv.Line((x_cog_f, 0, z_min), (x_cog_f, 0, z_max))
    line3 = pv.Line((x_cog, 0, z_min), (x_cog, 0, z_max))

    # Add a text block next to the line
    text_position = (x_cog_w + 0.2, 0, z_max)  # Slightly offset in x-direction
    text_position2 = (x_cog_f + 0.2, 0, z_max)
    text_position3 = (x_cog + 0.2, 0, z_max)
    if aircraft_type == "conventional":
        x_tail = x_cog + lv
        y_tail = 0
        z_tail = 2 * r_fus
        htail1 = plot_straight_wing(file_path_0012, span=ref.b_h/2, chord=ref.c_root_h, dihedral_deg=0,
                                    start_point=(x_tail + 0.5 *(ref.c_root_v-ref.c_tip_v), y_tail, z_tail + ref.b_v))
        vtail = plot_vertical_wing_t(file_path_0012, span=ref.b_v, root_chord=ref.c_root_v, tip_chord=ref.c_tip_v,
                                     sweep_deg=0, start_point=(x_tail, y_tail, z_tail))
        htail2 = plot_straight_wing(file_path_0012, span=-ref.b_h/2, chord=ref.c_root_h, dihedral_deg=0,
                                    start_point=(x_tail + 0.5 * (ref.c_root_v-ref.c_tip_v), y_tail, z_tail + ref.b_v))
        empennage = htail1 + vtail + htail2
        plotter.add_text("Reference Aircraft View", position='upper_edge', font_size=12, color='black')
    elif aircraft_type == "DUUC":
        duuc_left = visualize_propulsive_empennage_ac(pe_input[0], pe_input[1], pe_input[2], pe_input[3],
                                                      pe_input[4], pe_input[5], pe_input[6], pe_input[7],
                                                      pe_input[8], pe_input[9], pe_input[10], pe_input[11],
                                                      pe_input[12], (x_duuc, y_duuc, z_duuc))

        duuc_right = duuc_left.reflect((0, 1, 0))
        empennage = duuc_left + duuc_right
        plotter.add_text("DUUC Aircraft View", position='upper_edge', font_size=12, color='black')
    else:
        raise ValueError("Invalid aircraft type specified. Choose 'conventional' or 'DUUC'.")
    # Plot the result
    #plotter = pv.Plotter()
    plotter.show_axes()
    plotter.set_background("skyblue")
    plotter.add_mesh(body, color=off_white, show_edges=False, smooth_shading=True)
    plotter.add_mesh(wing_comb, color=off_white, show_edges=False, smooth_shading=True)
    plotter.add_mesh(empennage, color='lightgrey', show_edges=False, smooth_shading=True)

    # add center of gravity functionality
    plotter.add_point_labels([text_position], ['x_(cg-wing)'], font_size=14, point_color='white', point_size=10)
    plotter.add_mesh(line, color="black", line_width=2)
    plotter.add_point_labels([text_position2], ['x_(cg-fus)'], font_size=14, point_color='white', point_size=10)
    plotter.add_mesh(line2, color="black", line_width=2)
    plotter.add_point_labels([text_position3], ['x_(cg)'], font_size=14, point_color='white', point_size=10)
    plotter.add_mesh(line3, color="red", line_width=2)
    plotter.camera_position = [(-7, 0, 10), (2.5, 2.5, -0.5), (0, 0, 1)]
    plotter.camera.zoom(-12)
    plotter.show()


def screen_param():
    from screeninfo import get_monitors

    monitors = get_monitors()
    second_monitor = monitors[2]  # Second screen (index 1)

    # Get the position of the second monitor
    x_position = second_monitor.x
    y_position = second_monitor.y
    return x_position, y_position

"""
PE_input = [config.pylon_chord, config.pylon_length, config.duct_chord, config.duct_diameter,
            config.support_chord, config.support_length, config.cant_angle,
            config.control_vane_chord, config.control_vane_length, 0.25 * config.duct_chord,
            0.95 * config.duct_chord, 0.40 * config.duct_chord, 0.30 * config.duct_chord]

visualize_aircraft(3.7, ref.b_w/2, ref.c_root_w, 11.25, 11.5, 25, 15,
                   aircraft_type="DUUC", lv=14, pe_input=PE_input) """
