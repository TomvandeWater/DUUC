import matlab.engine
import numpy as np
from analysis_modules.ISA import air_density_isa
import flow_conditions, config
BEM_matlab_engine = matlab.engine.start_matlab()


def bem_sequence(n_blades, diameter, velocity, j_start, j_end, density, temperature, propname):
    BEM_matlab_engine.cd(r'C:\Users\tomva\pythonProject\DUUC\analysis_modules\BEM')
    n_iterations = int((j_end - j_start) * 10 + 1)

    advance = np.linspace(j_start, j_end, n_iterations)
    print(advance)
    BEM_vector = []
    for i in range(len(advance)):
        t_out, q_out, n_out, tc, cp, ct = BEM_matlab_engine.BEM2(n_blades, diameter, 10, 0, velocity, 0,
                                                                 advance[i], density, temperature, propname, nargout=6)
        BEM_vector.append([advance[i], t_out, q_out, n_out, tc, cp, ct])

    # Convert to NumPy array if needed
    BEM_vector = np.array(BEM_vector)

    return BEM_vector

"""
density_cruise = air_density_isa(flow_conditions.altitude)[0]
temperature_cruise = air_density_isa(flow_conditions.altitude)[2]
v_cruise = 128

thrustdata = bem_sequence(config.n_blades, config.duct_diameter, v_cruise, 0.5, 1.0, density_cruise, temperature_cruise, 'HM568F')
print(f"Thrust vector: {thrustdata}")"""
