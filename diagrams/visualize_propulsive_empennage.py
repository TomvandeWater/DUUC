import numpy as np
import os
import config
import pyvista as pv
pv.global_theme.allow_empty_mesh = True


logo_texture = pv.read_texture(r"C:\Users\tomva\pythonProject\DUUC\data\images\TU_Delft.png")


def cmyk_to_rgb(c, m, y, k):
    """
    Converts CMYK values to RGB values.

    Parameters:
        c (float): Cyan value (0 to 1).
        m (float): Magenta value (0 to 1).
        y (float): Yellow value (0 to 1).
        k (float): Key/Black value (0 to 1).

    Returns:
        tuple: RGB values (r, g, b), each between 0 and 255.
    """
    r = int(255 * (1 - c) * (1 - k))
    g = int(255 * (1 - m) * (1 - k))
    b = int(255 * (1 - y) * (1 - k))
    return (r, g, b)


def apply_single_logo_texture(mesh, image_path, center, size_u=0.5, size_v=0.5, direction_u=(0, 1, 0), direction_v=(0, 0, 1)):
    # Load texture image
    texture = pv.read_texture(image_path)
    texture.repeat = False  # Do not repeat across surface

    # Calculate U and V points for texture mapping
    point_u = np.array(center) + size_u * np.array(direction_u)
    point_v = np.array(center) + size_v * np.array(direction_v)

    # Apply plane mapping with controlled bounds
    textured_mesh = mesh.texture_map_to_plane(
        origin=center,
        point_u=point_u,
        point_v=point_v,
        inplace=False,
        use_bounds=False,
    )

    return textured_mesh, texture


def visualize_annular_wing(file_path, chord=1.0, diameter=2.0,
                           start_point=(0, 0, 0), num_sections=100):
    """
    Create an annular wing mesh from airfoil coordinates.

    Parameters:
    file_path (str): Path to airfoil coordinate file (2D X-Z points).
    chord (float): Scaling factor for airfoil profile.
    diameter (float): Diameter of the ring wing.
    start_point (tuple): Center point of the ring wing (x, y, z).
    num_sections (int): Number of cross-sections around the ring.

    Returns:
    pyvista.PolyData: Annular wing mesh.
    """
    # Load and prepare airfoil profile
    airfoil = load_airfoil(file_path, chord)

    # Rotate airfoil thickness by 90 degrees
    airfoil.points = rotate_thickness_90(airfoil.points)

    # Generate circular path parameters
    radius = diameter / 2
    theta = np.linspace(0, 2 * np.pi, num_sections)

    # Initialize points and faces for the annular wing mesh
    all_points = []
    all_faces = []

    for i, angle in enumerate(theta):
        # Calculate position on circular path (Y-Z plane)
        y_pos = start_point[1] + radius * np.cos(angle)
        z_pos = start_point[2] + radius * np.sin(angle)
        position = [start_point[0], y_pos, z_pos]

        # Calculate orientation vectors
        tangent = np.array([0, -np.sin(angle), np.cos(angle)])  # Tangent to circular path
        normal = np.array([1, 0, 0])  # Chordwise direction (X-axis)
        binormal = np.cross(tangent, normal)  # Perpendicular to tangent and normal

        # Create rotation matrix
        rot_matrix = np.column_stack((normal, binormal, tangent))  # Align X-Z plane

        # Rotate and translate airfoil
        section_points = airfoil.copy()
        section_points.points = section_points.points @ rot_matrix.T + position

        # Append points to the global list
        all_points.append(section_points.points)

        # Create faces connecting this section to the previous one
        if i > 0:
            prev_points = all_points[i - 1]
            curr_points = section_points.points
            for j in range(len(curr_points) - 1):
                all_faces.extend([
                    4,  # Quad face with 4 points
                    len(all_points[i - 1]) * (i - 1) + j,
                    len(all_points[i - 1]) * i + j,
                    len(all_points[i - 1]) * i + j + 1,
                    len(all_points[i - 1]) * (i - 1) + j + 1,
                ])

    # Flatten points array and create PolyData
    all_points_flattened = np.vstack(all_points)
    wing_mesh = pv.PolyData()
    wing_mesh.points = all_points_flattened
    wing_mesh.faces = np.array(all_faces)

    return wing_mesh


def load_airfoil(file_path, chord):
    """Load and prepare airfoil profile from file."""
    # Read 2D airfoil coordinates (X-Z plane)
    coords = np.loadtxt(file_path)
    x = coords[:, 0] * chord
    z = coords[:, 1] * chord

    # Ensure closed loop (connect last to first point)
    if not np.allclose(coords[0], coords[-1]):
        x = np.append(x, x[0])
        z = np.append(z, z[0])

    # Create 3D polyline (Y=0 initially; will be rotated later)
    points = np.c_[x, np.zeros_like(x), z]
    polyline = pv.PolyData()
    polyline.points = points
    polyline.lines = np.hstack([[len(points)], range(len(points))])

    return polyline


def rotate_thickness_90(points):
    """Rotate airfoil thickness by swapping Z and Y axes."""
    rotated_points = points.copy()
    rotated_points[:, [1, 2]] = rotated_points[:, [2, 1]]  # Swap Y and Z axes
    return rotated_points


def plot_straight_wing(file_path, span=1.0, chord=1.0, dihedral_deg=0.0, start_point=(0,0,0), color='lightblue', normalized=True):
    # Load airfoil coordinates
    airfoil_coords = np.loadtxt(file_path)
    x_raw, y_raw = airfoil_coords[:,0], airfoil_coords[:,1]

    if normalized:
        x2d = x_raw * chord
        y2d = y_raw * chord
    else:
        x2d = x_raw * chord
        y2d = y_raw  # Already in physical units

    # Airfoil in 3D: align X (chord), Y (span), Z (thickness)
    root_profile = np.column_stack((x2d, np.zeros_like(x2d), y2d))

    # Dihedral calculation
    half_span = span
    dihedral_rad = np.radians(dihedral_deg)
    tip_y_offset = half_span
    tip_z_offset = half_span * np.tan(dihedral_rad)

    # Tip profile (translated)
    tip_profile = root_profile + np.array([0, tip_y_offset, tip_z_offset])

    # Translate to start point
    root_profile += np.array(start_point)
    tip_profile += np.array(start_point)

    # Combine points
    all_points = np.vstack((root_profile, tip_profile))
    n_points = len(root_profile)

    # Quad faces between root and tip
    faces = []
    for i in range(n_points - 1):
        p1 = i
        p2 = i + n_points
        p3 = i + 1 + n_points
        p4 = i + 1
        faces.append([4, p1, p2, p3, p4])  # Quad

    # --- Closing the root section ---
    root_center = np.mean(root_profile, axis=0)
    root_center_idx = all_points.shape[0]
    all_points = np.vstack((all_points, root_center))
    for i in range(n_points - 1):
        p1 = i
        p2 = i + 1
        faces.append([3, root_center_idx, p1, p2])  # Triangle fan

    # --- Closing the tip section ---
    tip_center = np.mean(tip_profile, axis=0)
    tip_center_idx = all_points.shape[0]
    all_points = np.vstack((all_points, tip_center))
    for i in range(n_points - 1):
        p1 = i + n_points
        p2 = i + 1 + n_points
        faces.append([3, tip_center_idx, p2, p1])  # Triangle fan (reverse order)

    faces = np.hstack(faces)

    # Create mesh
    wing_mesh = pv.PolyData(all_points, faces)
    return wing_mesh


def plot_vertical_wing(file_path, span=1.0, chord=1.0, sweep_deg=0.0, start_point=(0,0,0), color='lightgreen', normalized=True):
    airfoil_coords = np.loadtxt(file_path)
    x_raw, y_raw = airfoil_coords[:,0], airfoil_coords[:,1]

    if normalized:
        x2d = x_raw * chord
        y2d = y_raw * chord
    else:
        x2d = x_raw * chord
        y2d = y_raw

    # Root profile: map airfoil Y to Y-axis; lies in X-Y, translated in Z
    root_profile = np.column_stack((x2d, y2d, np.zeros_like(x2d)))  # Z=0 at root

    # Sweep calculation: span in Z, offset in Y
    sweep_rad = np.radians(sweep_deg)
    tip_z_offset = span
    tip_y_offset = span * np.tan(sweep_rad)

    # Tip profile (translated in Z and Y)
    tip_profile = root_profile + np.array([0, tip_y_offset, tip_z_offset])

    # Translate to start point
    root_profile += np.array(start_point)
    tip_profile += np.array(start_point)

    # Combine points
    all_points = np.vstack((root_profile, tip_profile))
    n_points = len(root_profile)

    # Quad faces
    faces = []
    for i in range(n_points - 1):
        p1 = i
        p2 = i + n_points
        p3 = i + 1 + n_points
        p4 = i + 1
        faces.append([4, p1, p2, p3, p4])

    # Root cap
    root_center = np.mean(root_profile, axis=0)
    root_center_idx = len(all_points)
    all_points = np.vstack((all_points, root_center))
    for i in range(n_points - 1):
        p1 = i
        p2 = i + 1
        faces.append([3, root_center_idx, p1, p2])

    # Tip cap
    tip_center = np.mean(tip_profile, axis=0)
    tip_center_idx = len(all_points)
    all_points = np.vstack((all_points, tip_center))
    for i in range(n_points - 1):
        p1 = i + n_points
        p2 = i + 1 + n_points
        faces.append([3, tip_center_idx, p2, p1])

    faces = np.hstack(faces)

    wing_mesh = pv.PolyData(all_points, faces)

    return wing_mesh


def plot_half_sphere(diameter=1.0, center_point=(0, 0, 0)):
    # Create a full sphere with radius = diameter / 2
    sphere = pv.Sphere(radius=diameter / 2, center=(0, 0, 0))

    # Clip the sphere to get the half-sphere (keep the part where x >= 0)
    half_sphere = sphere.clip_box(bounds=(0, np.inf, -np.inf, np.inf, -np.inf, np.inf))

    half_sphere.translate(center_point, inplace=True)

    return half_sphere


def create_propeller(file_path, span, chord, num_wings, center_point=(0, 0, 0)):
    # Load airfoil coordinates
    airfoil_coords = np.loadtxt(file_path)
    x2d, z2d = airfoil_coords[:, 0], airfoil_coords[:, 1]

    # Normalize and scale to specified chord
    x2d = x2d - np.min(x2d)
    chord_length = np.max(x2d)
    if chord_length == 0:
        raise ValueError("Airfoil X coordinates have zero length.")
    x2d *= (chord / chord_length)
    z2d *= (chord / chord_length)

    # Create airfoil surface span-wise
    y_vals = np.linspace(0, span, 2)  # Two span-wise positions
    points = []
    faces = []

    for j, y in enumerate(y_vals):
        for i in range(len(x2d)):
            points.append([x2d[i], y, z2d[i]])

    n_points_per_section = len(x2d)
    for i in range(n_points_per_section - 1):
        faces.append([4, i, i + 1, i + 1 + n_points_per_section, i + n_points_per_section])

    points = np.array(points)
    faces = np.hstack(faces)

    airfoil_mesh = pv.PolyData(points, faces)

    # Create multiple rotated wings
    combined = pv.MultiBlock()
    angle_step = 360 / num_wings

    for n in range(num_wings):
        angle = n * angle_step
        rotated_mesh = airfoil_mesh.rotate_x(angle, point=(0, 0, 0), inplace=False)
        rotated_mesh.translate(center_point, inplace=True)
        combined.append(rotated_mesh)

    # Merge all wings into one mesh for output
    return combined.combine()


def create_engine_nacelle(diameter=1.0, length=2.0, cone_length=0.5, center_point=(0, 0, 0), resolution=100):
    radius = diameter / 2

    # Create the cylinder (aligned along X-axis)
    cylinder = pv.Cylinder(center=(length/2, 0, 0), direction=(1, 0, 0),
                           radius=radius, height=length, resolution=resolution)

    # Create the cone (aligned along X-axis)
    cone = pv.Cone(center=(length, 0, 0), direction=(1, 0, 0),
                   height=cone_length, radius=radius, resolution=resolution)

    # Combine the cone and cylinder
    nacelle = cylinder + cone

    # Translate so that the non-conical end starts at center_point
    nacelle.translate(center_point, inplace=True)

    return nacelle


def visualize_propulsive_empennage(plotter, c_pylon, b_pylon, c_duct, d_duct, c_support, b_support, cant, c_cv, b_cv, LE_pylon,
                                   LE_control, LE_support, x_prop, origin=(0, 0, 0)):

    primary = cmyk_to_rgb(1, 0, 0, 0)
    off_white = cmyk_to_rgb(0, 0, 0.02, 0.02)
    file_path_0012 = os.path.join(r"C:\Users\tomva\pythonProject\DUUC\data\airfoil_coordinates", "Naca0012" + ".txt")
    cant_rad = np.radians(cant)

    # Displace the origin by adding the origin offsets
    origin_x, origin_y, origin_z = origin

    # Create pylon, support, duct, and control surfaces with displacement
    pylon = plot_straight_wing(file_path_0012, span=np.cos(cant_rad) * b_pylon, chord=c_pylon, dihedral_deg=cant,
                               start_point=(LE_pylon + origin_x, origin_y, origin_z))

    x_pylon_end = LE_support + origin_x
    y_pylon_end = np.cos(cant_rad) * b_pylon + origin_y
    z_pylon_end = np.sin(cant_rad) * b_pylon + origin_z

    support = plot_straight_wing(file_path_0012, span=np.cos(cant_rad) * b_support, chord=c_support, dihedral_deg=cant,
                                 start_point=(x_pylon_end, y_pylon_end, z_pylon_end))

    x_duct = 0 + origin_x
    y_duct = y_pylon_end + np.cos(cant_rad) * (0.5 * d_duct) + origin_y
    z_duct = z_pylon_end + np.sin(cant_rad) * (0.5 * d_duct) + origin_z

    duct = visualize_annular_wing(file_path_0012, chord=c_duct, diameter=d_duct,
                                  start_point=(x_duct, y_duct, z_duct), num_sections=100)

    logo_width = c_duct * 0.75  # meters
    logo_height = c_duct * 0.4  # meters
    logo_center = (x_duct + 0.15 * c_duct, y_duct + d_duct / 2, z_duct - 0.5 * logo_height)  # Top of the ring

    # Apply logo texture to the duct mesh
    duct_with_logo, logo_texture = apply_single_logo_texture(
        duct,
        r"C:\Users\tomva\pythonProject\DUUC\data\images\TU_Delft.png",
        center=logo_center,
        size_u=logo_width,
        size_v=logo_height,
        direction_u=(1, 0, 0),  # Logo width runs along Z+
        direction_v=(0, 0, 1)  # Logo height runs along X+
    )

    x_control_h = LE_control + origin_x
    y_control_h = y_pylon_end + np.cos(cant_rad) * b_support / 2 - 0.5 * d_duct + origin_y
    z_control_h = z_pylon_end + np.sin(cant_rad) * b_support / 2 + origin_z
    h_control = plot_straight_wing(file_path_0012, span=b_cv * 2, chord=c_cv, dihedral_deg=0,
                                   start_point=(x_control_h, y_control_h, z_control_h))

    x_control_v = LE_control + origin_x
    y_control_v = y_pylon_end + np.cos(cant_rad) * b_support / 2 + origin_y
    z_control_v = z_pylon_end + np.sin(cant_rad) * b_support / 2 - 0.5 * d_duct + origin_z
    v_control = plot_vertical_wing(file_path_0012, span=b_cv * 2, chord=c_cv, sweep_deg=0.0,
                                   start_point=(x_control_v, y_control_v, z_control_v))

    x_engine = x_prop + origin_x
    spinner = plot_half_sphere(config.hub_diameter, (x_engine, y_duct, z_duct))

    propeller = create_propeller(file_path_0012, span=d_duct/2, center_point=(x_engine - config.c_root, y_duct, z_duct),
                                 num_wings=config.n_blades, chord=config.c_root)

    nacelle = create_engine_nacelle(diameter=config.hub_diameter, length=config.nacelle_length,
                                    cone_length=0.5, center_point=(x_engine, y_duct, z_duct), resolution=100)

    # Create primary plotter with the full geometry and table
    plotter.set_background("skyblue")

    # Add meshes to plotter, ensuring that everything is displaced by the origin
    plotter.add_text("Propulsive Empennage", position='upper_edge', font_size=12, color='black')
    plotter.add_mesh(pylon, color=off_white, smooth_shading=True)
    plotter.add_mesh(support, color=off_white, smooth_shading=True)
    plotter.add_mesh(h_control, color="darkblue", smooth_shading=True)
    plotter.add_mesh(v_control, color="darkblue", smooth_shading=True)
    plotter.add_mesh(duct_with_logo, texture=logo_texture, smooth_shading=True)
    plotter.add_mesh(duct, color="lightgrey", smooth_shading=True)
    plotter.add_mesh(spinner, color="grey", smooth_shading=True)
    plotter.add_mesh(propeller, color="grey", smooth_shading=True)
    plotter.add_mesh(nacelle, color="lightgrey", smooth_shading=True)
    plotter.camera_position = [(-7, 0, 10), (2.5, 2.5, -0.5), (0, 0, 1)]
    plotter.camera.zoom(-12)
    plotter.add_axes()

    plotter.show()


def visualize_cross_section(plotter, c_pylon, b_pylon, c_duct, d_duct, c_support, b_support, cant, c_cv, b_cv, LE_pylon,
                            LE_control, LE_support, x_prop):

    primary = cmyk_to_rgb(1, 0, 0, 0)
    file_path_0012 = os.path.join(r"C:\Users\tomva\pythonProject\DUUC\data\airfoil_coordinates", "Naca0012" + ".txt")
    cant_rad = np.radians(cant)

    pylon = plot_straight_wing(file_path_0012, span=np.cos(cant_rad) * b_pylon, chord=c_pylon, dihedral_deg=cant,
                               start_point=(LE_pylon, 0, 0))
    x_pylon_end = LE_support
    y_pylon_end = np.cos(cant_rad) * b_pylon
    z_pylon_end = np.sin(cant_rad) * b_pylon

    support = plot_straight_wing(file_path_0012, span=np.cos(cant_rad) * b_support, chord=c_support, dihedral_deg=cant,
                                 start_point=(x_pylon_end, y_pylon_end, z_pylon_end))

    x_duct = 0
    y_duct = y_pylon_end + np.cos(cant_rad) * (0.5 * d_duct)
    z_duct = z_pylon_end + np.sin(cant_rad) * (0.5 * d_duct)

    duct = visualize_annular_wing(file_path_0012, chord=c_duct, diameter=d_duct,
                                  start_point=(x_duct, y_duct, z_duct), num_sections=100)

    logo_width = c_duct * 0.75  # meters
    logo_height = c_duct * 0.4  # meters
    logo_center = (x_duct + 0.15 * c_duct, y_duct + d_duct / 2, z_duct - 0.5 * logo_height)  # Top of the ring

    # Apply logo texture to the duct mesh
    duct_with_logo, logo_texture = apply_single_logo_texture(
        duct,
        r"C:\Users\tomva\pythonProject\DUUC\data\images\TU_Delft.png",
        center=logo_center,
        size_u=logo_width,
        size_v=logo_height,
        direction_u=(1, 0, 0),  # Logo width runs along Z+
        direction_v=(0, 0, 1)  # Logo height runs along X+
    )

    x_control_h = LE_control
    y_control_h = y_pylon_end + np.cos(cant_rad) * b_support / 2 - 0.5 * d_duct
    z_control_h = z_pylon_end + np.sin(cant_rad) * b_support / 2
    h_control = plot_straight_wing(file_path_0012, span=b_cv * 2, chord=c_cv, dihedral_deg=0,
                                   start_point=(x_control_h, y_control_h, z_control_h))
    x_control_v = LE_control
    y_control_v = y_pylon_end + np.cos(cant_rad) * b_support / 2
    z_control_v = z_pylon_end + np.sin(cant_rad) * b_support / 2 - 0.5 * d_duct
    v_control = plot_vertical_wing(file_path_0012, span=b_cv * 2, chord=c_cv, sweep_deg=0.0,
                                   start_point=(x_control_v, y_control_v, z_control_v))

    x_engine = x_prop
    spinner = plot_half_sphere(config.hub_diameter, (x_engine, y_duct, z_duct))

    propeller = create_propeller(file_path_0012, span=d_duct / 2,
                                 center_point=(x_engine - config.c_root, y_duct, z_duct),
                                 num_wings=config.n_blades, chord=config.c_root)

    nacelle = create_engine_nacelle(diameter=config.hub_diameter, length=config.nacelle_length,
                                    cone_length=0.5, center_point=(x_engine, y_duct, z_duct), resolution=100)

    # --- Create secondary plotter for cross-sectional view with slider ---
    #cross_section_plotter = pv.Plotter(window_size=(800, 800), title="Cross-Sectional View")
    plotter.set_background("skyblue")
    #cross_section_plotter.add_text("Cross Section: XY Plane", position='upper_edge', font_size=12, color='black')

    # Add the original full geometry for reference (light grey, transparent)
    combined = pylon + support + h_control + v_control + duct_with_logo + spinner + propeller + nacelle
    plotter.add_mesh(combined, color="lightgrey", opacity=0.2)

    # Initial Z position for slicing
    initial_z = z_duct  # center of the duct location

    # Slice function to update with slider
    def update_slice(z_value):
        plotter.clear_actors()  # Clear previous slice and re-add reference
        plotter.add_mesh(combined, color="lightgrey", opacity=0.2)
        plotter.add_text("Cross Section: XY Plane", position='upper_edge', font_size=12, color='black')

        # Create the slice through all components at z=z_value
        slice_plane_origin = (0, 0, z_value)
        slice_normal = (0, 0, 1)  # Normal in Z-direction (XY-plane slice)

        sliced = combined.slice(origin=slice_plane_origin, normal=slice_normal)
        plotter.add_mesh(sliced, color="black", line_width=2)
        plotter.add_axes()
        plotter.view_xy()

    # Add slider to control Z slicing position
    z_min = z_duct - d_duct * 0.5  # Extend below duct
    z_max = z_duct + d_duct * 0.5  # Extend above duct
    plotter.add_slider_widget(update_slice, value=initial_z,
                                            rng=[z_min, z_max],
                                            title='Z-slider [m]',
                                            pointa=(0.225, .1), pointb=(0.425, .1),
                                            style='modern')

    # Initial slice display
    update_slice(initial_z)
    plotter.show()


def get_parameters(plotter, c_pylon, b_pylon, c_duct, d_duct, c_support, b_support, cant, c_cv, b_cv, LE_pylon,
                  LE_control, LE_support, x_prop):
    # parameters for table
    params = {
        "Pylon Chord": f"{c_pylon} [m]",
        "Pylon Span": f"{b_pylon} [m]",
        "Cant angle": f"{cant} [deg]",
        "Support Chord": f"{c_support} [m]",
        "Support Span": f"{b_support} [m]",
        "Duct Chord": f"{c_duct} [m]",
        "Duct Diameter": f"{d_duct} [m]",
        "Duct Aspect Ratio": f"{d_duct / c_duct} [-]",
        "Number of blades": f"{config.n_blades} [-]",
        "Leading Edge Duct": f"{0} [m]",
        "Leading Edge Pylon": f"{LE_pylon} [m]",
        "Leading Edge Support": f"{LE_support} [m]",
        "Leading Edge Control": f"{LE_control} [m]",
        "Propeller location": f"{x_prop} [m]",
    }

    plotter.add_text("Parameters", font_size=12, position="upper_edge")

    y_pos = 0.85
    for key, value in params.items():
        plotter.add_text(f"{key}: {value}", position=(0.05, y_pos), font_size=10, viewport=True)
        y_pos -= 0.04  # Space between lines

    plotter.show()


"""
visualize_cross_section(plotter, config.pylon_chord, config.pylon_length, config.duct_chord, config.duct_diameter,
                               config.support_chord, config.support_length, config.cant_angle,
                               config.control_vane_chord, config.control_vane_length, 0.25 * config.duct_chord,
                               0.95 * config.duct_chord, 0.40 * config.duct_chord, 0.30 * config.duct_chord)

plotter = pv.Plotter()
visualize_propulsive_empennage(plotter, config.pylon_chord, config.pylon_length, config.duct_chord, config.duct_diameter,
                               config.support_chord, config.support_length, config.cant_angle,
                               config.control_vane_chord, config.control_vane_length, 0.25 * config.duct_chord,
                               0.95 * config.duct_chord, 0.40 * config.duct_chord, 0.30 * config.duct_chord)"""

