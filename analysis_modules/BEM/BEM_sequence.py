import matlab.engine
import numpy as np
from analysis_modules.ISA import air_density_isa
import flow_conditions, config
import matplotlib.pyplot as plt
BEM_matlab_engine = matlab.engine.start_matlab()


def bem_sequence(n_blades, diameter, velocity, j_start, j_end, density, temperature, propname):
    BEM_matlab_engine.cd(r'C:\Users\tomva\pythonProject\DUUC\analysis_modules\BEM')
    n_iterations = int((j_end - j_start) * 10 + 1)

    advance = np.linspace(j_start, j_end, n_iterations)
    print(advance)
    BEM_vector = []
    for i in range(len(advance)):
        t_out, q_out, n_out, tc, cp, ct, Va_avg, Vt_avg, cn = BEM_matlab_engine.BEM2(n_blades, diameter, 36, 0, velocity, 5,
                                                                 advance[i], density, temperature, propname, nargout=9)
        BEM_vector.append([advance[i], t_out / 1000, q_out / 1000, n_out / 1000, tc, cp, ct, Va_avg, Vt_avg, cn])

    # Convert to NumPy array if needed
    BEM_vector = np.array(BEM_vector)

    return BEM_vector


density_cruise = air_density_isa(7000)[0]
temperature_cruise = -1
v_cruise = 50

output_vector = bem_sequence(config.n_blades, config.duct_diameter, v_cruise, 0.6, 1.6, density_cruise,
                             temperature_cruise, 'prop1')


def interpolate_coefficients(output_vector, requested_advance):
    advance_ratios = output_vector[:, 0]
    t_out_values = output_vector[:, 1]
    q_out_values = output_vector[:, 2]
    n_out_values = output_vector[:, 3]
    tc_values = output_vector[:, 4]
    cp_values = output_vector[:, 5]
    ct_values = output_vector[:, 6]
    Va_avg = output_vector[:, 7]
    Vt_avg = output_vector[:, 8]
    cn = output_vector[:, 9]

    # Perform linear interpolation using np.interp
    interpolated_t_out = np.interp(requested_advance, advance_ratios, t_out_values)
    interpolated_q_out = np.interp(requested_advance, advance_ratios, q_out_values)
    interpolated_n_out = np.interp(requested_advance, advance_ratios, n_out_values)
    interpolated_tc = np.interp(requested_advance, advance_ratios, tc_values)
    interpolated_cp = np.interp(requested_advance, advance_ratios, cp_values)
    interpolated_ct = np.interp(requested_advance, advance_ratios, ct_values)
    interpolated_va = np.interp(requested_advance, advance_ratios, Va_avg)
    interpolated_vt = np.interp(requested_advance, advance_ratios, Vt_avg)
    interpolated_cn = np.interp(requested_advance, advance_ratios, cn)


    return {
        'advance_ratio': requested_advance,
        't_out': interpolated_t_out,
        'q_out': interpolated_q_out,
        'n_out': interpolated_n_out,
        'tc': interpolated_tc,
        'cp': interpolated_cp,
        'ct': interpolated_ct
    }





advance = output_vector[:, 0]
t_out = output_vector[:, 1]
q_out = output_vector[:, 2]
n_out = output_vector[:, 3]
tc = output_vector[:, 4]
cp = output_vector[:, 5]
ct = output_vector[:, 6]
Va_avg = output_vector[:, 7]
Vt_avg = output_vector[:, 8]
cn = output_vector[:, 9]

# Plot
plt.figure("Forces vs. advance ratio", figsize=(10, 6))
plt.plot(advance, t_out, label='Thrust')
plt.plot(advance, q_out, label='Torque')
plt.plot(advance, n_out, label='Normal force')
plt.xlabel('Advance Ratio [-]')
plt.ylabel('Forces [kN]')
plt.title('Forces vs. Advance Ratio')
plt.legend()
plt.grid()

plt.figure("Coefficients vs. advance ratio", figsize=(10, 6))
plt.plot(advance, tc, label=r'$t_c$')
plt.plot(advance, cp, label=r'$c_p$')
plt.plot(advance, ct, label=r'$c_t$')
plt.plot(advance, cn, label=r'$c_n$')
plt.xlabel('Advance Ratio [-]')
plt.ylabel('Coefficients [-]')
plt.title('Coefficients vs. Advance Ratio')
plt.legend()
plt.grid()

plt.figure("Velocities vs. advance ratio", figsize=(10, 6))
plt.plot(advance, Va_avg, label=r'$V_{axial}$')
plt.plot(advance, Vt_avg, label=r'$V_{tangential}$')
plt.xlabel('Advance Ratio [-]')
plt.ylabel('Velocity [m/s]')
plt.title('Velocities vs. Advance Ratio')
plt.legend()
plt.grid()

plt.show()

