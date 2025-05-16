import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter

""" plot for control vane deflection effect on effective duct inflow angle"""

"""
aspect_ratio = 2
"""
filepath1 = r"C:\Users\tomva\pythonProject\DUUC\data\AVL\AVL_coeff_comp.txt"
filepath2 = r"C:\Users\tomva\pythonProject\DUUC\data\AVL\AVL_elev_deflect.txt"
"""

column_names = ['alpha', 'cl_duct', 'cd_duct', 'cm_duct', 'cl_cv1', 'cd_cv1', 'cm_cv1',
                'cl_support', 'cd_support', 'cm_support', 'cl_cv2', 'cd_cv2', 'cm_cv2']
df = pd.read_csv(filepath1, sep=r'\s+', skiprows=2, names=column_names)
column_names2 = ['deflect', 'cl_duct', 'cl_cv', 'cd_duct', 'cd_cv', 'cl_tot', 'cd_tot', 'cm_tot']
df2 = pd.read_csv(filepath2, sep=r'\s+', skiprows=2, names=column_names2)

alpha1 = df["cl_duct"]

cd1 = df["cd_duct"]
cl_de = df2["cl_duct"]  #  op c = 1.8
cd_de = df2["cd_duct"]
x_vect = np.linspace(-1, 1, 101)
de_2 = [0.0296, -0.0082, -0.0404, -0.0759]  # op c = 1.5

a = (1.3218 + 0.0082) / 3
b = -0.0082


# Define function for reuse
def alpha_from_deflection(delta_e_input):
    return slope * delta_e_input + intercept


def y(x):
    y_line = a * x + b
    return y_line


def scale_cl_labels(value, tick_number):
    return f"{value / 5:.2f}"


y_vect = [y(x) for x in x_vect]

# Plot the main line
plt.figure()
plt.plot(x_vect, y_vect, color="tab:blue", label='Duct lift coefficient')

colors = ['tab:orange', 'tab:red', 'tab:green', 'tab:purple']
labels = [r'$\delta_e$ = -5', r'$\delta_e$ = 0', r'$\delta_e$ = 5', r'$\delta_e$ = 10']

for i, cl in enumerate(cl_de):
    x_intersect = (cl - b) / a
    y_intersect = cl
    plt.hlines(y=cl, xmin=-1, xmax=x_intersect, color=colors[i], linestyle='dashed', label=labels[i])
    plt.plot(x_intersect, y_intersect, 'o', color=colors[i])
    plt.vlines(x_intersect, ymin=0, ymax=y_intersect, linestyle='dashed', color=colors[i])
    delta_e = df2["deflect"].to_numpy()
    alpha_from_intersect = [(cl - b) / a for cl in cl_de]
    alpha_from_intersect = np.array(alpha_from_intersect)
    slope, intercept = np.polyfit(delta_e, alpha_from_intersect, 1)
    print(f"Empirical relation: alpha = {slope:.4f} * delta_e + {intercept:.4f}")

# Plot x-markers for each y-value in de_2 on the blue line
for i, y_val in enumerate(de_2):
    x_val = (y_val - b) / a
    plt.plot(x_val, y_val, 'x', color=colors[i], markersize=8, markeredgewidth=2)


plt.ylim(-0.15, 0.1)
plt.xlim(-0.5, 0.5)
plt.gca().yaxis.set_major_formatter(FuncFormatter(scale_cl_labels))
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$C_L$')
plt.title('Elevator Deflection Effect on Duct Inflow Angle')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()
plt.show()"""

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import data.experiment_ct_duct as ct_duct
from matplotlib.transforms import Bbox
import data.experiment_reference_7foot as ref7

# Load data
advance1 = ct_duct.advance_naca0012
ct1 = ct_duct.ct_duct_naca0012
cd_a0 = 0.011363

cd_duct_emp = [0.05258286047262917, 0.02828497097369604, 0.02447822164241713, 0.022521593224637322, 0.021248476932469385, 0.02032304196198435, 0.01960509707022396, 0.019023667722862594, 0.018538201563751446, 0.018123506647080248, 0.01776293966144966, 0.017444972575364762, 0.017161315225591448, 0.016905821200279385, 0.016673816661130027, 0.016461670856118243, 0.01626651138389712, 0.01608602964717858, 0.015918344429443898, 0.01576190404009626, 0.015615414713711162, 0.015477787287459709, 0.015348096860748851, 0.015225551841919741, 0.01510946989254111, 0.014999259014368521, 0.01489440252163009, 0.01479444698440809, 0.014698992469330867, 0.014607684574803357, 0.014520207881303895, 0.014436280527296326, 0.01435564968779959, 0.014278087782304768, 0.014203389276173227, 0.014131367968155171, 0.01406185467856057, 0.013994695269565996, 0.013929748942369511, 0.013866886766303347, 0.01380599040324315, 0.013746950997207684, 0.013689668204296339, 0.01363404934234671, 0.013580008643127398, 0.013527466592677994, 0.013476349347698717, 0.013426588217776588, 0.013378119204792534, 0.013330882592146525, 0.013284822577515431, 0.013239886943759892, 0.013196026763353833, 0.013153196132348626, 0.013111351930424016, 0.013070453604036185, 0.013030462970063673, 0.01299134403768515, 0.012953062846508618, 0.01291558731921681, 0.012878887127204978, 0.012842933567869638, 0.012807699452365054, 0.012773159002781176, 0.012739287757816383, 0.01270606248612232, 0.01267346110658922, 0.012641462614919772, 0.012610047015909565, 0.012579195260913823, 0.012548889190034251, 0.012519111478607752, 0.012489845587621264, 0.012461075717714407, 0.012432786766465085, 0.01240496428868267, 0.012377594459460094, 0.012350664039759481, 0.012324160344327099, 0.012298071211752422, 0.012272384976502652, 0.012247090442779586, 0.012222176860059053, 0.01219763390018558, 0.012173451635905936, 0.012149620520735205, 0.012126131370057995, 0.012102975343375682, 0.012080143927617776, 0.012057628921442408, 0.01203542242045682, 0.012013516803294428, 0.01199190471848996, 0.011970579072098775, 0.011949533016010675, 0.011928759936912276, 0.011908253445855565, 0.011888007368393348, 0.011868015735245346, 0.011848272773461241, 0.011828772898049488, 0.011809510704042982, 0.011790480958974635, 0.011771678595737948, 0.011753098705809333, 0.011734736532810514, 0.011716587466390973, 0.011698647036411564, 0.011680910907411877, 0.01166337487334501, 0.011646034852564437, 0.011628886883048818, 0.011611927117851343, 0.011595151820761197, 0.011578557362165445, 0.011562140215100417, 0.011545896951482377]
advance_emp = [0.0015503875968992248, 0.01704089815557338, 0.03253140871424753, 0.04802191927292168, 0.06351242983159583, 0.07900294039026998, 0.09449345094894414, 0.10998396150761829, 0.12547447206629242, 0.1409649826249666, 0.15645549318364074, 0.1719460037423149, 0.18743651430098907, 0.2029270248596632, 0.21841753541833736, 0.2339080459770115, 0.24939855653568568, 0.26488906709435983, 0.28037957765303395, 0.2958700882117081, 0.31136059877038225, 0.3268511093290564, 0.34234161988773054, 0.3578321304464047, 0.3733226410050789, 0.38881315156375307, 0.4043036621224272, 0.41979417268110136, 0.4352846832397755, 0.45077519379844966, 0.4662657043571238, 0.48175621491579795, 0.49724672547447213, 0.5127372360331462, 0.5282277465918204, 0.5437182571504946, 0.5592087677091687, 0.5746992782678428, 0.590189788826517, 0.6056802993851912, 0.6211708099438653, 0.6366613205025394, 0.6521518310612135, 0.6676423416198878, 0.6831328521785619, 0.698623362737236, 0.7141138732959101, 0.7296043838545844, 0.7450948944132586, 0.7605854049719327, 0.7760759155306068, 0.7915664260892811, 0.8070569366479552, 0.8225474472066293, 0.8380379577653034, 0.8535284683239777, 0.8690189788826518, 0.8845094894413259, 0.9, 0.9154905105586743, 0.9309810211173484, 0.9464715316760225, 0.9619620422346966, 0.9774525527933708, 0.9929430633520449, 1.008433573910719, 1.023924084469393, 1.0394145950280673, 1.0549051055867413, 1.0703956161454156, 1.0858861267040898, 1.1013766372627638, 1.116867147821438, 1.1323576583801123, 1.1478481689387863, 1.1633386794974605, 1.1788291900561347, 1.1943197006148087, 1.209810211173483, 1.225300721732157, 1.2407912322908312, 1.2562817428495054, 1.2717722534081795, 1.2872627639668537, 1.302753274525528, 1.318243785084202, 1.3337342956428762, 1.3492248062015502, 1.3647153167602244, 1.3802058273188986, 1.3956963378775726, 1.4111868484362469, 1.426677358994921, 1.442167869553595, 1.4576583801122693, 1.4731488906709433, 1.4886394012296178, 1.504129911788292, 1.519620422346966, 1.5351109329056403, 1.5506014434643145, 1.5660919540229885, 1.5815824645816627, 1.5970729751403367, 1.612563485699011, 1.6280539962576852, 1.6435445068163592, 1.6590350173750334, 1.6745255279337077, 1.6900160384923817, 1.705506549051056, 1.72099705960973, 1.7364875701684042, 1.7519780807270784, 1.7674685912857524, 1.7829591018444266, 1.7984496124031009]


advance2 = ct_duct.advance_naca0018
ct2 = ct_duct.ct_duct_naca0018

advance3 = ct_duct.advance_naca4312
ct3 = ct_duct.ct_duct_naca4312

# Polynomial fits
poly1 = np.polyfit(advance1, ct1, deg=2)
poly2 = np.polyfit(advance2, ct2, deg=2)
poly3 = np.polyfit(advance3, ct3, deg=2)

p1 = np.poly1d(poly1)
p2 = np.poly1d(poly2)
p3 = np.poly1d(poly3)
print(f"p1: {p1}")
print(f"p2: {p2}")
print(f"p3: {p3}")
# Adjust fit ranges to start at J=0.35
start_advance = 0.35  # Starting advance ratio for extrapolation
advance_fit1 = np.linspace(start_advance, 1.7, 100)
advance_fit2 = np.linspace(start_advance, 1.7, 100)
advance_fit3 = np.linspace(start_advance, 1.7, 100)

# Main plot
fig, ax = plt.subplots()
line1, = ax.plot(advance1, ct1, label='Naca0012', marker='o')
line2, = ax.plot(advance2, ct2, label='Naca0018', marker='o')
line3, = ax.plot(advance3, ct3, label='Naca4312', marker='o')


# Dashed extrapolated fits starting at J=0.35
ax.plot(advance_fit1, p1(advance_fit1), linestyle='--', color=line1.get_color())
ax.plot(advance_fit2, p2(advance_fit2), linestyle='--', color=line2.get_color())
ax.plot(advance_fit3, p3(advance_fit3), linestyle='--', color=line3.get_color())
ax.plot(advance_emp, [-cd for cd in cd_duct_emp], label='Emperical prediction')

# Extrapolated points at 1.7
ax.plot(1.7, p1(1.7), color=line1.get_color())
ax.plot(1.7, p2(1.7), color=line2.get_color())
ax.plot(1.7, p3(1.7), color=line3.get_color())
#ax.plot([0, 1.75], [-cd_a0, -cd_a0], label='Duct drag - empirical', linestyle='dashed')

# Labels and grid
ax.set_title('Thrust coefficient of the Duct')
ax.set_ylabel(r'$C_{T,duct}$ [-]')
ax.set_xlabel('Advance Ratio [-]')
ax.set_ylim(top=0.05)
ax.grid(True)
ax.legend(loc='upper right')

# Inset axes: zoomed view up to J=0.4, offset to avoid axis labels
bbox = Bbox.from_bounds(-0.35, -0.50, 1, 1)  # Adjusted position and size
axins = inset_axes(ax, width="40%", height="40%", bbox_to_anchor=bbox, bbox_transform=ax.transAxes)

# Plot same data in inset
axins.plot(advance1, ct1, marker='o', color=line1.get_color())
axins.plot(advance2, ct2, marker='o', color=line2.get_color())
axins.plot(advance3, ct3, marker='o', color=line3.get_color())

# Inset limits
axins.set_xlim(0, 0.4)
axins.set_ylim(bottom=min(min(ct1), min(ct2), min(ct3)) * 0.95, top=0.011)
axins.grid(True)

# Draw connector lines
mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")


ct_duct_interp = np.interp(ct_duct.advance_prop_naca4312, ct_duct.advance_naca4312, ct_duct.ct_duct_naca4312)

# Now the arrays have same length and can be added
ct_total = ct_duct.ct_prop_naca4312 + ct_duct_interp

# Create the plot
plt.figure()
plt.plot(ref7.ct_j, ref7.ct, label='7 foot ct')
plt.plot(ct_duct.advance_prop_naca4312, ct_total, label='CT prop naca4312')
plt.plot(advance_emp, [-cd for cd in cd_duct_emp], label='emperical')


# Add grid, legend, and show plot
plt.grid(True)
plt.title("Thrust Coefficient Ducted Fan")
plt.ylabel(r"$C_T$ [-]")
plt.xlabel(r"$J [-]$")
plt.legend()
plt.show() """


elev_deflect = [-10, -5, 0, 5, 10]
elev_box = [-4.38, -2.19, 0, 2.19, 4.38]
elev_rw = [-1.1881, -0.5941, 0, 0.5941, 1.1881]

elev_box_conv = []
elev_rw_conv = []
for i in range(len(elev_deflect)):
    elev_box_conv.append(elev_box[i] / 3.5)
    elev_rw_conv.append(elev_rw[i] / 3.5)

grad_box = (elev_box_conv[4] - elev_box_conv[0]) / 20
grad_rw = (elev_rw_conv[4] - elev_rw_conv[0]) / 20

print(f"grad box: {grad_box}, grad rw: {grad_rw}")

grad_box_rad = grad_box * 180 / np.pi
grad_rw_rad = grad_rw * 180 / np.pi

plt.figure()
plt.plot(elev_deflect, [elev_box / 3.5 for elev_box in elev_box], label='Elevator on Duct Edge')
plt.plot(elev_deflect, [elev_rw / 3.5 for elev_rw in elev_rw], label='Elevator V0.1 Configuration')

plt.title('Elevator Effectiveness')
plt.ylabel(r"$C_N$ [-]")
plt.xlabel(r"$\delta_e$ [deg]")
plt.legend()
plt.grid(True)
plt.show()

