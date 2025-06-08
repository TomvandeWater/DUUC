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

"""
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



grad_box_rad = grad_box * 180 / np.pi
grad_rw_rad = grad_rw * 180 / np.pi

print(f"grad box: {grad_box_rad}, grad rw: {grad_rw_rad}")

plt.figure()
plt.plot(elev_deflect, [elev_box / 3.5 for elev_box in elev_box], label='Elevator on Duct Edge')
plt.plot(elev_deflect, [elev_rw / 3.5 for elev_rw in elev_rw], label='Elevator V0.1 Configuration')

plt.title('Elevator Effectiveness')
plt.ylabel(r"$C_N$ [-]")
plt.xlabel(r"$\delta_e$ [deg]")
plt.legend()
plt.grid(True)


import data.experiment_reference_5annular_airfoil as ar

plt.figure()
plt.plot(ar.cla_a_ar_1_3, ar.cla_cl_ar_1_3, label='AR = 1/3', marker='o')
plt.plot(ar.cla_a_ar_2_3, ar.cla_cl_ar_2_3, label='AR = 2/3', marker='s')
plt.plot(ar.cla_a_ar_1_0, ar.cla_cl_ar_1_0, label='AR = 1.0', marker='D')
plt.plot(ar.cla_a_ar_1_5, ar.cla_cl_ar_1_5, label='AR = 1.5', marker=">")
plt.plot(ar.cla_a_ar_3_0, ar.cla_cl_ar_3_0, label='AR = 3.0', marker='^')

plt.title(r'$C_L$ vs. $\alpha$ for different Aspect Ratios')
plt.ylabel(r"$C_L$ [-]")
plt.xlabel(r"$\alpha$ [deg]")
plt.legend()
plt.grid(True)
plt.show()"""
"""
import data.experiment_reference_prop_location as prop

a = np.linspace(0, 15, 16)
cma = []
AR =1.62

xp = 0.3 * AR ** -0.768 + - 0.26
xe = 0.63 * AR + -0.35
kp = 6.25 * np.sin(AR / 2)
kv = np.pi / 3

for i in range(len(a)):
    alpha = np.radians(a[i])
    cma.append(xp * kp * np.sin(alpha) * np.cos(alpha) + xe * kv * np.sin(alpha) ** 2)

plt.figure()
plt.plot(a, cma, label='Theor.', color='black')
plt.scatter(prop.maq_cma_a_sc, prop.maq_cma_cm_sc, label="Exp.", color='black')
plt.ylabel(r"$C_M$ [-]")
plt.xlabel(r"$\alpha$ [deg]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()"""

"""
cd_atr = [0.027405602923264313, 0.025578562728380026, 0.02466504263093788, 0.02484774665042631, 0.025213154689403168, 0.026492082825822168, 0.027953714981729598, 0.029780755176613886, 0.03233861144945189, 0.03526187576126675, 0.03818514007308161, 0.0412911084043849, 0.044762484774665046, 0.04841656516443362, 0.05225334957369062, 0.05627283800243606, 0.060109622411693056, 0.06467722289890378, 0.06887941534713764, 0.0730816077953715, 0.07764920828258222, 0.08239951278928137, 0.08678440925700366]
cl_atr = [-0.1991701244813278, -0.08713692946058091, 0.006224066390041493, 0.09647302904564314, 0.1991701244813278, 0.29253112033195017, 0.3890041493775933, 0.4730290456431535, 0.5539419087136929, 0.6410788381742738, 0.7219917012448133, 0.7935684647302904, 0.8651452282157676, 0.9398340248962654, 1.008298755186722, 1.070539419087137, 1.1327800829875518, 1.2012448132780082, 1.2572614107883817, 1.3070539419087137, 1.3692946058091284, 1.4221991701244812, 1.475103734439834]
cd_duuc1 = [0.025213154689403168, 0.023934226552984165, 0.023568818514007307, 0.023751522533495738, 0.025761266747868453, 0.028501827040194886, 0.03361753958587089, 0.03818514007308161, 0.04275274056029233, 0.04823386114494519, 0.05261875761266748, 0.05682095006090134, 0.0610231425091352, 0.0654080389768575, 0.0712545676004872, 0.07710109622411693, 0.0825822168087698, 0.09190012180267966, 0.10249695493300853]
cl_duuc1 = [-0.2520746887966805, -0.16182572614107882, -0.08402489626556016, 0.07468879668049792, 0.23651452282157676, 0.37655601659751037, 0.529045643153527, 0.6473029045643153, 0.7406639004149377, 0.8371369294605808, 0.9180497925311203, 0.9740663900414938, 1.0331950207468878, 1.0954356846473028, 1.1701244813278007, 1.2417012448132778, 1.300829875518672, 1.3973029045643153, 1.5]
cd_duuc2 = [0.025943970767356883, 0.024299634591961022, 0.023934226552984165, 0.023751522533495738, 0.025030450669914737, 0.026857490864799025, 0.028501827040194886, 0.03087697929354446, 0.033252131546894034, 0.03635809987819732, 0.03928136419001218, 0.04275274056029233, 0.046041412911084045, 0.04987819732034105, 0.05353227771010963, 0.05791717417783192, 0.06175395858708892, 0.06613885505481121, 0.07034104750304507, 0.07490864799025579, 0.07965895249695494, 0.08440925700365408, 0.08989037758830695, 0.09482338611449452]
cl_duuc2 = [-0.1929460580912863, -0.09024896265560166, 0.012448132780082987, 0.11203319502074688, 0.20850622406639002, 0.2987551867219917, 0.38278008298755184, 0.4730290456431535, 0.5539419087136929, 0.6317427385892116, 0.7095435684647302, 0.7748962655601659, 0.8464730290456431, 0.9056016597510372, 0.9678423236514522, 1.0394190871369293, 1.0985477178423235, 1.154564315352697, 1.2074688796680497, 1.266597510373444, 1.3226141078838174, 1.3692946058091284, 1.437759336099585, 1.496887966804979]

cd_duuc1 = [0.025580046403712295, 0.024883990719257537, 0.024187935034802783, 0.023491879350348025, 0.023317865429234336, 0.023491879350348025, 0.023491879350348025, 0.024013921113689093, 0.024361948955916472, 0.024883990719257537, 0.025754060324825984, 0.026624129930394428, 0.027320185614849186, 0.02853828306264501, 0.029930394431554524, 0.031148491879350346, 0.03236658932714617, 0.033584686774941995, 0.034976798143851504, 0.036194895591647326, 0.037412993039443156, 0.038805104408352664, 0.04054524361948956, 0.04211136890951276, 0.04367749419953596, 0.04524361948955916, 0.04698375870069605, 0.048375870069605566, 0.05011600928074245, 0.05237819025522041, 0.05464037122969837, 0.05655452436194895, 0.05846867749419953, 0.06038283062645011, 0.06229698375870069, 0.06455916473317865, 0.0668213457076566, 0.06925754060324825, 0.07134570765661252, 0.07360788863109048, 0.0765661252900232, 0.07987238979118329, 0.08387470997679813, 0.08735498839907192, 0.08996519721577725, 0.09205336426914153, 0.0941415313225058, 0.09622969837587006, 0.0986658932714617, 0.1004060324825986]
cl_duuc1 = [-0.2589285714285714, -0.21726190476190474, -0.15773809523809523, -0.09226190476190475, -0.047619047619047616, -0.002976190476190476, 0.03869047619047619, 0.09226190476190475, 0.14285714285714285, 0.1875, 0.22916666666666666, 0.2648809523809524, 0.30952380952380953, 0.3630952380952381, 0.40773809523809523, 0.4464285714285714, 0.4851190476190476, 0.5178571428571428, 0.5595238095238095, 0.5892857142857143, 0.6190476190476191, 0.6517857142857143, 0.6875, 0.7232142857142857, 0.7529761904761905, 0.7827380952380952, 0.8154761904761905, 0.8422619047619048, 0.869047619047619, 0.9107142857142857, 0.9494047619047619, 0.9761904761904762, 1.005952380952381, 1.0297619047619047, 1.0595238095238095, 1.0922619047619047, 1.1160714285714286, 1.1488095238095237, 1.169642857142857, 1.2053571428571428, 1.2380952380952381, 1.2708333333333333, 1.3154761904761905, 1.3541666666666665, 1.3779761904761905, 1.4017857142857142, 1.425595238095238, 1.4464285714285714, 1.4732142857142856, 1.4970238095238095]
cd_atr = [0.027320185614849186, 0.026450116009280742, 0.025754060324825984, 0.025232018561484916, 0.025058004640371227, 0.024535962877030162, 0.024535962877030162, 0.024883990719257537, 0.024883990719257537, 0.025058004640371227, 0.025754060324825984, 0.026798143851508117, 0.027668213457076565, 0.02853828306264501, 0.029234338747099766, 0.030800464037122968, 0.03149651972157772, 0.03306264501160092, 0.034628770301624125, 0.03602088167053364, 0.03706496519721578, 0.03897911832946636, 0.040371229698375866, 0.042285382830626446, 0.04385150812064965, 0.04576566125290023, 0.04733178654292343, 0.0494199535962877, 0.050986078886310904, 0.053596287703016235, 0.05498839907192575, 0.057424593967517396, 0.05916473317865429, 0.061600928074245935, 0.06334106728538283, 0.06577726218097447, 0.06751740139211136, 0.06995359628770301, 0.07204176334106728, 0.07430394431554524, 0.07674013921113688, 0.07900232018561484]
cl_atr = [-0.20238095238095238, -0.15178571428571427, -0.11011904761904762, -0.05952380952380952, -0.011904761904761904, 0.03869047619047619, 0.0744047619047619, 0.13392857142857142, 0.17857142857142855, 0.22023809523809523, 0.2708333333333333, 0.31845238095238093, 0.3571428571428571, 0.4107142857142857, 0.44940476190476186, 0.5, 0.5297619047619048, 0.5773809523809523, 0.6190476190476191, 0.6666666666666666, 0.6994047619047619, 0.7410714285714285, 0.7738095238095237, 0.8184523809523809, 0.8511904761904762, 0.8898809523809523, 0.9196428571428571, 0.9583333333333333, 0.9910714285714285, 1.0327380952380951, 1.056547619047619, 1.0922619047619047, 1.1220238095238095, 1.1577380952380951, 1.1785714285714286, 1.2172619047619047, 1.2380952380952381, 1.2708333333333333, 1.294642857142857, 1.3273809523809523, 1.357142857142857, 1.380952380952381]

from scipy.signal import savgol_filter

def smooth_savgol(x, y, window_length=5, polyorder=2):
    x = np.array(x)
    y = np.array(y)
    if len(x) < window_length:
        window_length = len(x) if len(x) % 2 == 1 else len(x)-1
    y_smooth = savgol_filter(y, window_length=window_length, polyorder=polyorder)
    return x, y_smooth
from scipy.interpolate import CubicSpline
t = np.arange(len(cd_atr))  # or use cumulative distance for arc-length-like parameterization

# Create spline as function of t
spline_cd = CubicSpline(t, cd_atr)
spline_cl = CubicSpline(t, cl_atr)

# Evaluate
t_new = np.linspace(t[0], t[-1], 300)
cd_new = spline_cd(t_new)
cl_new = spline_cl(t_new)


# Smooth each curve
cd_atr_s, cl_atr_s = smooth_savgol(cd_atr, cl_atr)
cd_duuc1_s, cl_duuc1_s = smooth_savgol(cd_duuc1, cl_duuc1)
cd_duuc2_s, cl_duuc2_s = smooth_savgol(cd_duuc2, cl_duuc2)


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

x = np.array([0.03156976744186046, 0.035755813953488376, 0.03959302325581395, 0.042906976744186046, 0.04709302325581395, 0.05162790697674419, 0.05843023255813953])
y = np.array([0.47023809523809523, 0.5833333333333333, 0.6726190476190476, 0.7410714285714285, 0.8214285714285714, 0.8988095238095237, 1.0029761904761905])


atr_x = np.array([0.029825581395348837, 0.03244186046511628, 0.03488372093023256, 0.03784883720930233, 0.041511627906976746, 0.044651162790697675, 0.04848837209302326, 0.0518])
atr_y = np.array([0.46428571428571425, 0.5625, 0.6339285714285714, 0.7142857142857142, 0.7976190476190476, 0.869047619047619, 0.9375, 1.0])

# Generate new x values for a smooth curve
x_smooth = np.linspace(x.min(), x.max(), 500)
x_atr_smooth = np.linspace(atr_x.min(), atr_x.max(), 500)

# Create spline function
spline = make_interp_spline(x, y, k=3)  # k=3 cubic spline
spline_atr = make_interp_spline(atr_x, atr_y, k=3)  # k=3 cubic spline
# Generate smooth y values
y_smooth = spline(x_smooth)
y_smooth_atr = spline(x_atr_smooth)

def get_smooth_function(x, y):
    x = np.array(x)
    y = np.array(y)
    spline_function = make_interp_spline(x, y, k=3)
    return spline_function


# Plot
plt.figure()
plt.grid(True)
plt.plot(atr_x, atr_y, label="ATR 72-600", color="tab:orange")
plt.plot(x, y, label="DUUC - rear", color="tab:blue")
plt.plot(cd_duuc2_s, cl_duuc2_s, label="DUUC - front", color="tab:green")
plt.plot(0.03693, 0.69, color='orange', linestyle='--', label=r"$C_{L,cruise}$ ATR72-600", marker='o')
plt.plot(0.04139, 0.71, color='tab:blue', linestyle='--', label=r"$C_{L,cruise}$ DUUC - rear", marker='o')
plt.plot(0.03821, 0.68, color='tab:green', linestyle='--', label=r"$C_{L,cruise}$ DUUC - front", marker='o')
plt.title("Drag Polars")
plt.ylabel(r"$C_L$ [-]")
plt.xlabel(r"$C_D$ [-]")
plt.ylim(0.5, 1)
plt.legend()
"""

# Corrected Tail volume coefficients (scaled down)
CVv_values = [0.06, 0.07, 0.08, 0.09, 0.10, 0.11]
linestyles = ['-', '--', ':', '-', '--', ':']

# x-axis: l_h / b
lv_b = np.linspace(0.2, 0.6, 100)

# Initialize plot
plt.figure(figsize=(8, 6))

# Plot each line for different CVh
for CVv, style in zip(CVv_values, linestyles):
    Sh_S = CVv / lv_b
    plt.plot(lv_b, Sh_S, label=f'$C_{{V_h}}$ = {CVv:.2f}', linestyle=style)

# Add blue horizontal lines showing typical Sh/S range
plt.axhline(0.2, color='blue', linestyle=':', linewidth=1)
plt.axhline(0.36, color='blue', linestyle=':', linewidth=1)


# Labels and legend
plt.xlabel(r'$l_v/b$ [-]')
plt.ylabel(r'$S_v/S_w$ [-]')
plt.plot(0.5, 0.227, label="ATR72-600", marker='o', color='tab:orange')
plt.plot(0.35, 0.201, label="DUUC - front", marker='o', color='tab:blue')
plt.plot(0.255, 0.212, label="DUUC - rear", marker='o', color='tab:green')
plt.title("Vertical Tail Volume Coefficients")
plt.legend()
plt.grid(True)
plt.ylim(0.1, 0.4)
plt.xlim(0.2, 0.6)

# Corrected Tail volume coefficients (scaled down)
CVh_values = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
linestyles = ['-', '--', ':', '-', '--', ':']

# x-axis: l_h / b
lh_b = np.linspace(0.2, 7, 100)

# Initialize plot
plt.figure(figsize=(8, 6))

# Plot each line for different CVh
for CVh, style in zip(CVh_values, linestyles):
    Sh_S = CVh / lh_b
    plt.plot(lh_b, Sh_S, label=f'$C_{{V_h}}$ = {CVh:.2f}', linestyle=style)

# Add blue horizontal lines showing typical Sh/S range
plt.axhline(0.15, color='blue', linestyle=':', linewidth=1)
plt.axhline(0.30, color='blue', linestyle=':', linewidth=1)


# Labels and legend
plt.xlabel(r'$l_h/c$ [-]')
plt.ylabel(r'$S_h/S_w$ [-]')
plt.plot(6.04, 0.17, label="ATR72-600", marker='o', color='tab:orange')
plt.plot(4.25, 0.21, label="DUUC - front", marker='o', color='tab:blue')
plt.plot(3.08, 0.226, label="DUUC - rear", marker='o', color='tab:green')
plt.title("Horizontal Tail Volume Coefficients")
plt.legend()
plt.ylim(0.05, 0.40)
plt.xlim(2, 7)
plt.grid(True)
plt.show()

a_exp = [0, 1.0049019607843137, 1.9852941176470589, 3.002450980392157, 3.982843137254902, 5, 6.0171568627450975,
         6.985294117647059, 8.002450980392156, 9.019607843137255, 10.012254901960784, 10.992647058823529,
         11.985294117647058, 12.990196078431373, 13.995098039215685, 15]
cm_exp = [0.0019090909090909091, -0.005181818181818182, -0.009000000000000001, -0.01240909090909091,
          -0.015136363636363637, -0.01622727272727273, -0.01690909090909091, -0.020454545454545454,
          -0.02209090909090909, -0.024136363636363636, -0.024272727272727272, -0.02209090909090909,
          -0.020181818181818183, -0.018954545454545457, -0.01731818181818182, -0.014181818181818183]

a = np.linspace(0, 15, 32)
def cm(a):
    alfa = np.radians(a)
    ac = 0.3
    b = -0.786
    c = - 0.26
    AR = 1.62
    p1 = 0.63
    p0 = -0.35

    kpcm = 6.25 * np.sin(AR / 2)
    kvcm = np.pi / 3

    xp = ac * AR ** b + c
    xe = p1 * AR + p0

    cm_duct = xp * kpcm * np.sin(alfa) * np.cos(alfa) + xe * kvcm * np.sin(alfa) ** 2
    return cm_duct

plt.figure()
plt.plot(a, [cm(i) for i in a], label='Prediction model')
plt.plot(a_exp, cm_exp, label='Experimental data', marker='o', linestyle='dashed', color='tab:green')

plt.title(r'$C_m$ vs. $\alpha$ - Duct')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_m$ [-]')
plt.legend()
plt.grid(True)

import matplotlib.pyplot as plt
import numpy as np

# Your data
x = [0.5762304921968787, 2.2184873949579833, 3.8895558223289317, 5.675870348139256,
     7.404561824729892, 8.009603841536615, 8.672268907563025, 9.853541416566626,
     10.516206482593038, 11.582232893157263, 12.677070828331333, 13.599039615846339,
     14.117647058823529, 14.63625450180072, 15.471788715486195, 15.932773109243698,
     16.33613445378151, 16.797118847539014, 17.286914765906364, 17.690276110444177,
     18.036014405762305, 18.410564225690276, 18.842737094837936, 18.986794717887154]

y = [0.00007233273056057865, 0.0015189873417721517, 0.004122965641952983, 0.011934900542495477,
     0.02003616636528029, 0.02379746835443038, 0.02871609403254973, 0.037106690777576855,
     0.0420253164556962, 0.05070524412296564, 0.06083182640144666, 0.06922242314647378,
     0.07616636528028933, 0.08108499095840868, 0.09150090415913201, 0.09902350813743219,
     0.10625678119349005, 0.11233273056057866, 0.11985533453887884, 0.12853526220614828,
     0.13634719710669077, 0.14415913200723326, 0.14994575045207956, 0.15573236889692588]

# Convert to NumPy arrays
x = np.array(x)
y = np.array(y)

coeffs = np.polyfit(x, y, deg=2)
quadratic = np.poly1d(coeffs)

# Generate smooth x values
x_fit = np.linspace(min(x), max(x), 500)
y_fit = quadratic(x_fit)

a_cm_nac = [1.0660264105642256, 2.0168067226890756, 3.0252100840336134, 4.062424969987995, 4.638655462184874,
            5.7046818727490995, 5.647058823529412, 7.779111644657863, 7.779111644657863, 9.680672268907562,
            9.709483793517407, 11.78391356542617, 11.841536614645857, 13.97358943577431, 13.97358943577431,
            15.183673469387754, 16.30732292917167, 16.24969987995198, 17.229291716686674, 19.53421368547419]
cm_cm_nac = [0.00036166365280289325, 0.004122965641952983, 0.004990958408679927, 0.007305605786618444,
             0.009041591320072331, 0.00962025316455696, 0.014538878842676312, 0.017721518987341773,
             0.028426763110307417, 0.028426763110307417, 0.039710669077757686, 0.05649186256781194,
             0.06227848101265823, 0.09265822784810127, 0.09555153707052441, 0.10249547920433996,
             0.13403254972875225, 0.13027124773960216, 0.1490777576853526, 0.1916094032549729]

plt.figure()
plt.plot(x_fit, y_fit, label='Prediction model')
plt.scatter(a_cm_nac, cm_cm_nac, label='Experimental data', marker='o', linestyle='dashed', color='tab:green')

plt.title(r'$C_m$ vs. $\alpha$ - Nacelle')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_m$ [-]')
plt.legend()
plt.grid(True)

a1 = - 0.4887
b1 = 0.20768
a2 = 0.305
cmac = 2.2345
xcg = 0.3086
static_margin = 0.05
s_w = 61


x1 = np.linspace(-1, 0, 101)
x2 = np.linspace(-1, 2, 101)
x_ac = 0.25 * cmac

x_cg_ac = x1 + x_ac
y1 = a1 * x_cg_ac * cmac + b1
y2 = a2 * x2
y3 = y2 + static_margin

y_intersection = a2 * xcg + static_margin
x_intersect = (y_intersection - b1) / (a1 * cmac) - x_ac

sh_sw = np.round(y_intersection, 4)
new_sh = y_intersection * s_w
cg_range = xcg - x_intersect


a1_mod = - 0.4723
b1_mod = 0.23768
a2_mod = 0.29
cmac_mod = 2.2345
xcg_mod = 0.3286
static_margin_mod = 0.05
s_w = 61


x1_mod = np.linspace(-1, 0, 101)
x2_mod = np.linspace(-1, 2, 101)
x_ac_mod = 0.25 * cmac

x_cg_ac_mod = x1_mod + x_ac_mod
y1_mod = a1_mod * x_cg_ac_mod * cmac_mod + b1_mod
y2_mod = a2_mod * x2_mod
y3_mod = y2_mod + static_margin_mod

y_intersection_mod = a2_mod * xcg_mod + static_margin_mod
x_intersect_mod = (y_intersection_mod - b1_mod) / (a1_mod * cmac_mod) - x_ac_mod

sh_sw_mod = np.round(y_intersection_mod, 4)
new_sh_mod = y_intersection_mod * s_w
cg_range_mod = xcg_mod - x_intersect_mod


plt.figure(figsize=(12.8, 4.8))
plt.plot(x1, y1, label="Control - reference", linestyle='dashed')
plt.plot(x2, y2, label='Stability - reference', linestyle='dashed')
plt.plot(x2, y3, label='SM - reference', linestyle='dashed')

plt.plot([xcg, xcg], [0, y_intersection], linestyle="dotted", color='black', marker='>', label='cg - reference')
plt.plot([x_intersect, xcg], [y_intersection, y_intersection], linestyle="dotted", color='black')

plt.plot(x1_mod, y1_mod, label="Control - prediction model", color='tab:blue')
plt.plot(x2_mod, y2_mod, label='Stability - prediction model', color='tab:orange')
plt.plot(x2_mod, y3_mod, label='SM - prediction model', color='tab:green')
plt.plot([xcg_mod, xcg_mod], [0, y_intersection_mod], linestyle="dotted", color='purple', marker='>', label='cg - prediction model')
plt.plot([x_intersect_mod, xcg_mod], [y_intersection_mod, y_intersection_mod], linestyle="dotted", color='purple')


plt.title("Horizontal tail sizing - validation (zoomed)")
plt.ylabel(r"$S_h/S_w$ [-]")
plt.xlabel("Center of Gravity")
plt.xlim(-0.55, 0.55)
plt.ylim(0.1, 0.2)
plt.grid(True)
plt.legend()

array = np.linspace(0.001, 3, 301)


def s(a, b, c, d):
    s_final = 61 * (a/b) * (c/d)
    return s_final



plt.figure('CY_beta - vs S_V - validation')
plt.plot(array, [s(0.1914, t, 27.32, 13.27) for t in array], label=r'Prediction model', color="tab:blue")
plt.plot(array, [s(0.1811, t, 27.32, 14.43) for t in array], label=r'Reference', color="tab:blue", linestyle='dashed')

plt.plot([0, 2.228], [10.7, 10.7], color="tab:purple", linestyle="dashed")
plt.plot(2.228, 10.7, 'o', markersize=6, color="tab:purple", label="Prediction model value")
plt.plot([2.228, 2.228], [0, 10.7], linestyle='dashed', color="tab:purple")
plt.plot([0, 2.228], [9.57, 9.57], color="tab:orange", linestyle="dashed")
plt.plot(2.228, 9.57, 'o', markersize=6, color="tab:orange", label="Reference value")
plt.plot([2.228, 2.228], [0, 9.57], linestyle='dashed', color="tab:orange")
plt.xlim(0, 3)
plt.ylim(0, 40)
plt.xlabel(r'$C_{Y_{\beta}}$ [-]')
plt.ylabel(r'$S_{V}$ [$m^2$]')
plt.title(r'$S_{V}$ vs. $C_{Y_{\beta}}$ - validation')
plt.legend()
plt.grid(True)


def s_v(b, c):
    s_final = (51601 * np.pi * 2) / (0.5 * 1.225 * 65.43**2 * c * b)
    return s_final

plt.figure('CY - vs S_V - validation')
plt.plot(array, [s_v(t, 9.13) for t in array], label=r'Reference', color="tab:blue", linestyle='dashed')
plt.plot(array, [s_v(t, 8.75) for t in array], label=r'Prediction model', color="tab:blue", )
plt.plot([0, 0.925], [15.2, 15.2], color="tab:purple", linestyle="dashed")
plt.plot(0.925, 15.2, 'o', markersize=6, color="tab:purple", label="Prediction model value")
plt.plot([0.925, 0.925], [0, 15.2], linestyle='dashed', color="tab:purple")
plt.plot([0, 0.975], [14.085, 14.085], color="tab:orange", linestyle="dashed")
plt.plot(0.975, 14.085, 'o', markersize=6, color="tab:orange", label="Reference value")
plt.plot([0.975, 0.975], [0, 14.085], linestyle='dashed', color="tab:orange")



plt.xlim(0, 2)
plt.ylim(0, 40)
plt.xlabel(r'$C_{Y}$ [-]')
plt.ylabel(r'$S_{V}$ [$m^2$]')
plt.title(r'$S_{V}$ vs. $C_{Y}$ - validation')
plt.legend()
plt.grid(True)


a = [0, 5, 7, 10]
cl_1 = [0.6410, 1.2370, 1.4267, 1.4737]
cl_2 = [0, 0.6150, 0.8424, 1.0573]

cl_diff = []

for i in range(len(a)):
    cl_diff.append(cl_1[i] - cl_2[i])

plt.figure()
plt.plot(a, cl_diff, label="Correction difference", color='red', linestyle='dashed')
plt.plot(a, cl_1, label="Uncorrected control vane model")
plt.plot(a, cl_2, label="Corrected control vane model")

plt.xlabel(r'Elevator deflection $\delta_e$ [deg]')
plt.ylabel(r'$C_l$ [-]')
plt.title(r'Elevator lift coefficient comparison')
plt.legend()
plt.grid(True)



a_vik_poff = [-5, 5, 15, 20]
cl_vik_poff = [-0.3953033268101761, 0.42661448140900193, 1.0724070450097847, 0.9941291585127201]

a_dfs = np.linspace(0, 15, 31)
cl_dfs = [0.02775804748918083, 0.032487019003097715, 0.03721600025592129, 0.04194499492364383,
          0.04667197769950077, 0.051398880736767925, 0.056123779235729, 0.06084860534957858,
          0.06557278693359293,
          0.07029697048902829,
          0.07502389850407913,
          0.07975100335069503,
          0.08447680203722667,
          .08920258388607741,
          0.09391489301656138,
          0.09862658960680344,
          0.10339177725940056,
          0.10815967642275691,
          0.11304600143693844,
           0.11793826149630714,
          0.12254383440055291,
          .1271352818181481,
          0.13172681234675665,
                    0.13631842963203106,
                    0.14026559697361615,
                    0.14418093071987448,
                    0.14767974424466732,
                    0.15115802169113188,
                    0.1545133124121754,
                    0.15786261427043716,
                    0.1611159897479597]
cd_dfs = [0.01767399688743328,
                    0.01775820694149053,
                    0.01804810059896258,
                    0.018542566804443526,
                    0.019242800322584033,
                    0.020147617694097168,
                    0.021258133542824107,
                    0.022573240459815996,
                    0.0240935120912601,
                    0.025818362440601657,
                    0.027748628342622008,
                    0.029883502941128463,
                    0.03222381328629603,
                    0.034768754341961444,
                    0.037519545455317826,
                    0.04047501227257666,
                    0.04363705947537793,
                    0.047003846481194655,
                    0.050579197358707154,
                    0.05435941765982611,
                    0.058352687999867106,
                    0.06255108479425672,
                    0.0669542137300476,
                    0.07156209532554918,
                    0.07640934735040657,
                    0.0814631095179468,
                    0.0867951792554163,
                    0.09233573297533387,
                    0.09814917668741771,
                    0.1041708844878997,
                    0.110409774706254]

plt.figure()
plt.plot(a_dfs, [cl_dfs * 7 for cl_dfs in cl_dfs], label="Prediction model")
plt.plot(a_vik_poff, cl_vik_poff, label="Experimental", color='tab:green', marker='o')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_l$ [-]')
plt.title(r'Lift curve DFS - power off')
plt.legend()
plt.grid(True)

a_vik_pon = [-5, -0.023696682464454888, 5, 15, 20]
cl_vik_pon = [-0.8597194388777555, -0.14428857715430862, 0.721442885771543, 2.0861723446893787, 2.404809619238477]

cla_a_tc138_rpm2590 = [0.5, 9, 20]
cla_cl_tc138_rpm2590 = [0.7257019438444925, 3.83585313174946, 7.982721382289417]

plt.figure()
plt.plot(a_dfs, [cl_dfs * 14 -0.2 for cl_dfs in cl_dfs], label="Prediction model")
plt.plot(a_vik_pon, cl_vik_pon, label="Experimental Harinarain", color='tab:green', marker='o')
plt.plot(cla_a_tc138_rpm2590, [cl / 8.2 for cl in cla_a_tc138_rpm2590], label="Experimental Mort", marker='o')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_l$ [-]')
plt.title(r'Lift curve DFS - power on (J=0.42)')
plt.legend()
plt.grid(True)


clcd_cd_tc138_rpm2590 = [-13.846153846153847, -13.333333333333334, -12.41025641025641, -10.35897435897436, -8,
                         -4.102564102564104, -0.41025641025641235, 1.8461538461538467]
clcd_cl_tc138_rpm2590 = [0, 3.732181425485961, 8.086393088552915, 12.233261339092872, 15.44708423326134,
                         19.1792656587473, 21.771058315334773, 22.496760259179265]

plt.figure()
plt.plot([cd * 14 for cd in cd_dfs],
         [cl_dfs * 14 for cl_dfs in cl_dfs], label="Prediction model")
plt.plot([(cd / 10) + 1.555 for cd in clcd_cd_tc138_rpm2590],
         [cl / 10 for cl in clcd_cl_tc138_rpm2590], label="Experimental", color='tab:green', marker='o')
plt.xlabel(r'$C_d$ [-]')
plt.ylabel(r'$C_l$ [-]')
plt.title(r'Drag polar DFS - power-on ($T_c$=1.38)')
plt.legend()
plt.grid(True)

a_adib = [2.0567375886524824, 0.14184397163120568, 3.936170212765958, 5.99290780141844, 8.01418439716312, 10.070921985815604, 11.914893617021278, 14.007092198581562, 15.957446808510639, 18.049645390070925, 20]
cm_adib = [-0.0013483146067415732, -0.00006741573033707866, -0.002449438202247191, -0.0036404494382022475, -0.00453932584269663, -0.004921348314606742, -0.005303370786516854, -0.005235955056179776, -0.004853932584269663, -0.003910112359550562, -0.002359550561797753]

cm_dfs =[0.10761953578668623,
                    0.10662263179799637,
                    0.10562085107361807,
                    0.10461361577810668,
                    0.103604720193026,
                    0.10259093939413098,
                    0.1015755771798115,
                    0.1005557455007696,
                    0.09953406041170296,
                    0.0985082989167792,
                    0.09747803841388511,
                    0.09644393747714698,
                    0.09541050417363953,
                    0.09437390274383217,
                    0.0933523150123784,
                    0.09232880754734438,
                    0.09124853814169732,
                    0.09016262126671185,
                    0.08895125492765275,
                    0.0877300476239173,
                    0.08682366703308445,
                    0.08593239069569576,
                    0.08504070274651977,
                    0.08414879357788516,
                    0.08395429213256102,
                    0.08379965801782978,
                    0.0840634152937903,
                    0.084355383106604,
                    0.08476981260130859,
                    0.08519732886938722,
                    0.08574121426344211]

plt.figure()
plt.plot(a_dfs,
         [(cl_dfs * 0.8) - 0.08 for cl_dfs in cm_dfs], label="Prediction model")
plt.plot(a_adib, cm_adib, label="Experimental", color='tab:green', marker='o')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_m$ [-]')
plt.title(r'Pitching moment coefficient - validation')
plt.legend()
plt.grid(True)
import data.atr_reference as ref

a_atr = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15]
cl_atr = [0, 0.1248, 0.2493, 0.3731, 0.4966, 0.6199, 0.7421, 0.8649, 1.0088, 1.1592, 1.3316, 1.2733, 1.1950, 1.1039]
aratio = ref.s_ht / 61

a_pe = np.linspace(0, 15, 31)


plt.figure()
plt.plot(a_pe, [cl * 2 for cl in cl_dfs], label="Propulsive Empennage - power-off", color="tab:blue")
plt.plot(a_pe, [(cl * 1.18) * 2 for cl in cl_dfs], label="Propulsive Empennage - power-on", linestyle='dashed', color="tab:blue")
plt.plot(a_atr, [cl * aratio for cl in cl_atr], label="ATR72-600 empennage", color='tab:orange')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_l$ [-]')
plt.title(r'Lift curve - empennage comparison')
plt.legend()
plt.grid(True)

cd_atr = [0.00533, 0.00540, 0.00564, 0.00603, 0.00657, 0.00721, 0.00799, 0.00896, 0.01029, 0.01206, 0.02638, 0.04575, 0.07258, 0.10068]
e = 0.90


plt.figure()
#plt.plot([cd - 0.0118 for cd in cd_dfs], [cl * 2 for cl in cl_dfs], label="Propulsive Empennage - power-off", color="tab:blue")
plt.plot([cd * 1.18 - 0.0118*1.21 for cd in cd_dfs], [(cl * 1.18) * 2 for cl in cl_dfs], label="Propulsive Empennage" , color="tab:blue")
plt.plot([(0.8347e-3 + 1.315e-3) + ((cl ** 2) / (np.pi * e * ref.ar_h)) * aratio for cl in cl_atr], [cl * aratio for cl in cl_atr], label="ATR72-600 empennage", color='tab:orange')
plt.xlabel(r'$C_d$ [-]')
plt.ylabel(r'$C_l$ [-]')
plt.title(r'Drag polar - empennage comparison')
plt.legend()
plt.grid(True)

cm_atr = [0, 0.0004, 0.0009, 0.0018, 0.0033, 0.0054, 0.0085, 0.0122, 0.0124, 0.0124, 0.0498, 0.0434, 0.0269, 0.0115]
c_ratio = ref.c_mac_w / ref.c_root_h


plt.figure()
plt.plot(a_dfs, [cm * 2 for cm in cm_dfs], label="Propulsive Empennage - power-off", color="tab:blue")
plt.plot(a_dfs, [(cm * 1.18) * 2 for cm in cm_dfs], label="Propulsive Empennage - power-on", linestyle='dashed', color="tab:blue")
plt.plot(a_atr, [cm * aratio * c_ratio for cm in cm_atr], label="ATR72-600 empennage", color='tab:orange')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_m$ [-]')
plt.title(r'Pitching moment coefficient - empennage comparison')
plt.legend()
plt.grid(True)


class PropulsiveEmpennage:
    def __init__(self, l_pylon, cant):
        super().__init__()
        # operating conditions
        self.geometry_pylon = [l_pylon, 1, "0012", cant]
        self.geometry_support = [3.8, 1, "0012"]
        self.density = 0.548
        self.inflow_velocity = 128

    """ -------------------------------------- Weight ------------------------------------------------------------- """
    def weight_ps(self):
        # Safety factors
        n_ult = 1.5  # Ultimate load factor
        k_stoot = 1.5  # Impact factor
        k_misc = 1.25  # Miscellaneous factor

        # Convert cant angle to radians
        c_rad = np.radians(self.geometry_pylon[3])

        # Total length of the beam
        l_tot = self.geometry_support[0] + self.geometry_pylon[0]  # [m]

        # Material properties -> AL 6061-T6
        sigma_allow = 241 * 1e6  # Allowable stress [Pa] (converted from MPa to N/m^2)
        rho = 2700  # Mass per unit length of the beam [density]

        # Aerodynamic forces on pylon and support
        f_pylon = (0.5 * self.density * self.inflow_velocity ** 2 * 0.5 * self.geometry_pylon[1]
                   * self.geometry_pylon[0])  # [N]
        f_support = (0.5 * self.density * self.inflow_velocity ** 2 * 0.5
                     * self.geometry_support[1] * self.geometry_support[0])  # [N]

        # Aerodynamic moment about the root
        m_aero = (0.5 * f_pylon * self.geometry_pylon[0] + (0.5 * self.geometry_support[0] + self.geometry_pylon[0]) * f_support) * np.cos(
            c_rad)  # [N·m]

        # Weight moment about the root
        m_weight = (261+228+921) * (self.geometry_pylon[0] + self.geometry_support[0] * 0.5) * 9.81 * np.cos(
            c_rad)  # [N·m]

        # Sizing moment (maximum of aerodynamic and weight moments)
        m_sizing = max(m_aero, m_weight)  # [N·m]

        # Determine maximum thickness of the airfoil section
        num_list = [int(digit) for digit in self.geometry_pylon[2]]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness percentage
        thickness_pylon = (thickness / 100) * self.geometry_pylon[1]  # [m]

        num_list2 = [int(digit) for digit in self.geometry_support[2]]
        thickness_sup = num_list2[2] * 10 + num_list2[3]  # NACA thickness percentage
        thickness_support = (thickness_sup / 100) * self.geometry_support[1]  # [m]

        height = max(thickness_pylon, thickness_support) * 0.90  # Reduce maximum thickness by a safety factor [m]

        # Calculate width using bending stress formula
        width = (((6 * m_sizing * n_ult * k_stoot) / sigma_allow) / (height ** 2)) * k_misc  # [m]

        # Total mass of the beam
        mass_tot = height * width * l_tot * rho  # [kg]

        # Split mass into pylon and support sections
        #m_pylon = mass_of_section(mass_tot, l_tot, 0, self.geometry_pylon[0])  # Mass of pylon section [kg]
        #m_support = mass_of_section(mass_tot, l_tot, self.geometry_pylon[0], l_tot)  # Mass of support section [kg]

        return mass_tot

import data.atr_reference as ref
import config as config
aircraft_type ="DUUC"
sweep25 = ref.phi_qc_w
l_v = 9.5
power =2051 * 10 ** 3
eta = 0.9
approach_velocity = config.v_approach
n_engines = 2


def s_control(aircraft_type, sweep25, l_v, power, eta, approach_velocity, n_engines, y_engine, cy_duuc=0, cd_pe=0.00,
              cd_wind=0.00, cy_atr=0):
    velocity_mc = (1.3 / 1.2) * approach_velocity
    delta_f_max = 25
    swp = np.radians(sweep25)
    rho_sea = 1.225

    if aircraft_type == 'conventional':
        k_theta = (1 - 0.08 * np.cos(swp) ** 2) * np.cos(swp) ** (3/4)
        cl_the = 4.52  # from textNita -> still make this an interpollation based on input
        k_prime = 0.675
        cld_ratio = 0.85

        ne = ((eta * power) / (velocity_mc * n_engines)) * y_engine
        nd = 0.25 * ne
        top = ne + nd

        delta_f = np.radians(delta_f_max)
        denominator = (0.5 * rho_sea * velocity_mc ** 2 * delta_f * cld_ratio * cl_the * k_prime * k_theta * l_v)

    elif aircraft_type == 'DUUC':
        ne = ((eta * power) / (velocity_mc * n_engines)) * y_engine
        nd = (((cd_pe + cd_wind) * 0.5 * rho_sea * velocity_mc ** 2 * ref.s_w) * y_engine)
        top = ne + nd
        denominator = 0.5 * rho_sea * velocity_mc ** 2 * l_v * cy_duuc
    else:
        raise ValueError("Invalid aircraft type specified. Choose 'conventional' or 'DUUC'.")

    s_vertical = (top * 2 * np.pi) / denominator

    return s_vertical, top


y_engine = np.linspace(2.0, 5, 31)
moment = []
moment2 = []
for k in range(len(y_engine)):
    moment.append(s_control(aircraft_type, sweep25, l_v, power, eta, approach_velocity, n_engines, y_engine[k], cy_duuc=0.1,
                            cd_wind=0.06, cd_pe=0.026)[1])
    moment2.append(s_control(aircraft_type, sweep25, 6.5, power, eta, approach_velocity, n_engines, y_engine[k], cy_duuc=0.1,
                            cd_wind=0.06, cd_pe=0.022)[1])

plt.plot(np.linspace(0.2, 3, 31), moment, color="tab:blue")
plt.plot(1, 79604, label="DUUC V0.1", marker='o', color='tab:blue')
plt.xlabel(r'$Pylon\ length\ [m]$')
plt.ylabel(r'$Moment\ [Nm]$')
plt.title(r'Restoring moment required vs Pylon Length')
plt.legend()
plt.grid(True)
plt.tight_layout()

length = np.linspace(0.2, 3, 151)

def lift(length, cl):
    area = length * np.cos(np.radians(30)) * 1
    lift = cl * (2* area / 61)
    return lift


plt.figure()
plt.plot(length, [lift(length, 0.0) for length in length], color="tab:blue", label=r'$\alpha_{pylon}$ = 0')
plt.plot(length, [lift(length, 0.6251) for length in length], color="tab:orange", label=r'$\alpha_{pylon}$ = 5')
plt.plot(length, [lift(length, 1.1520) for length in length], color="tab:green", label=r'$\alpha_{pylon}$ = 10')
plt.axvline(x=1.0, label="DUUC V0.1", color="tab:blue", linestyle='dashed')
plt.xlabel(r'$Pylon\ length\ [m]$')
plt.ylabel(r'$C_L$ [-]')
plt.title(r'Lift contribution of the pylon vs Pylon Length')
plt.legend()
plt.grid(True)
plt.tight_layout()
"""
"""
def calculate_tail_area(
    x_cg_norm, x_tail_norm, fuselage_length, wing_area, mac, static_margin,
    x_ac_wing=0.25, eta_t=0.9, downwash_gradient=00,
    CL_alpha_w=5.406, CL_alpha_t=5.24, min_l_t=2.0
):
    # Convert normalized positions to actual positions (in meters)
    x_cg = x_cg_norm * fuselage_length
    x_tail = x_tail_norm * fuselage_length
    x_ac_wing_abs = x_ac_wing * mac

    # Tail moment arm
    l_t = x_tail - x_cg

    # Reject configurations where tail is unrealistically close to CG
    if abs(l_t) < min_l_t:
        return np.nan, l_t

    # Desired neutral point location
    x_np = x_cg + static_margin * mac

    # Tail contribution to neutral point (correct sign)
    tail_contribution = (x_np - x_ac_wing_abs) / (
        eta_t * (1 - downwash_gradient) * (CL_alpha_t / CL_alpha_w)
    )

    # Calculate required tail area, always positive
    S_t = abs(tail_contribution * wing_area * mac / l_t)
    s_t = np.sqrt(S_t)
    return s_t, l_t


x_tail = np.linspace(0.01, 1.0, 100)
x_cg = [9.2535791278755, 9.339134873554933, 9.424690619234369, 9.510246364913804, 9.59580211059324, 9.681357856272673, 9.766913601952108, 9.852469347631544, 9.938025093310976, 10.023580838990412, 10.109136584669848, 10.194692330349284, 10.280248076028716, 10.365803821708152, 10.451359567387588, 10.536915313067023, 10.622471058746456, 10.708026804425891, 10.793582550105327, 10.87913829578476, 10.964694041464195, 11.050249787143631, 11.135805532823067, 11.2213612785025, 11.306917024181935, 11.39247276986137, 11.478028515540807, 11.563584261220239, 11.649140006899675, 11.73469575257911, 11.820251498258543, 11.905807243937979, 11.991362989617414, 12.07691873529685, 12.162474480976286, 12.248030226655718, 12.333585972335154, 12.41914171801459, 12.504697463694022, 12.590253209373458, 12.675808955052894, 12.76136470073233, 12.846920446411762, 12.932476192091197, 13.018031937770633, 13.103587683450066, 13.189143429129501, 13.274699174808937, 13.360254920488373, 13.445810666167805, 13.531366411847241, 13.616922157526677, 13.702477903206113, 13.788033648885545, 13.87358939456498, 13.959145140244416, 14.044700885923849, 14.130256631603284, 14.21581237728272, 14.301368122962156, 14.386923868641588, 14.472479614321024, 14.55803536000046, 14.643591105679896, 14.729146851359328, 14.814702597038764, 14.9002583427182, 14.985814088397632, 15.071369834077068, 15.156925579756503, 15.24248132543594, 15.328037071115372, 15.413592816794807, 15.499148562474243, 15.584704308153679, 15.670260053833111, 15.755815799512547, 15.841371545191983, 15.926927290871415, 16.01248303655085, 16.098038782230287, 16.183594527909722, 16.269150273589155, 16.354706019268587, 16.440261764948026, 16.525817510627462, 16.611373256306894, 16.696929001986327, 16.782484747665766, 16.868040493345198, 16.95359623902463, 17.03915198470407, 17.124707730383506, 17.210263476062938, 17.29581922174237, 17.38137496742181, 17.466930713101245, 17.552486458780677, 17.638042204460113, 17.72359795013955]
x_cg2 = [7.2745525078163435, 7.341095865567015, 7.407639223317687, 7.474182581068358, 7.540725938819031, 7.607269296569702, 7.673812654320373, 7.740356012071044, 7.806899369821716, 7.873442727572387, 7.93998608532306, 8.006529443073731, 8.073072800824402, 8.139616158575073, 8.206159516325744, 8.272702874076415, 8.339246231827087, 8.405789589577761, 8.472332947328432, 8.538876305079103, 8.605419662829775, 8.671963020580446, 8.738506378331117, 8.805049736081788, 8.871593093832463, 8.938136451583134, 9.004679809333805, 9.071223167084476, 9.137766524835147, 9.204309882585818, 9.27085324033649, 9.33739659808716, 9.403939955837835, 9.470483313588506, 9.537026671339177, 9.603570029089848, 9.67011338684052, 9.73665674459119, 9.803200102341862, 9.869743460092533, 9.936286817843204, 10.002830175593875, 10.06937353334455, 10.13591689109522, 10.202460248845892, 10.269003606596563, 10.335546964347234, 10.402090322097905, 10.468633679848576, 10.535177037599247, 10.601720395349922, 10.668263753100593, 10.734807110851264, 10.801350468601935, 10.867893826352606, 10.934437184103281, 11.000980541853949, 11.067523899604623, 11.13406725735529, 11.200610615105965, 11.267153972856637, 11.333697330607308, 11.400240688357979, 11.46678404610865, 11.533327403859325, 11.599870761609996, 11.666414119360667, 11.732957477111338, 11.799500834862009, 11.86604419261268, 11.932587550363351, 11.999130908114022, 12.065674265864693, 12.132217623615364, 12.19876098136604, 12.26530433911671, 12.331847696867381, 12.398391054618052, 12.464934412368724, 12.531477770119398, 12.598021127870066, 12.66456448562074, 12.731107843371408, 12.797651201122083, 12.864194558872754, 12.930737916623425, 12.997281274374096, 13.063824632124767, 13.130367989875442, 13.19691134762611, 13.263454705376784, 13.329998063127455, 13.396541420878126, 13.463084778628797, 13.529628136379468, 13.59617149413014, 13.66271485188081, 13.729258209631485, 13.795801567382153, 13.862344925132827]

x_cg_norm = []
x_cg_norm2 =[]
for j in range(len(x_cg)):
    x_cg_norm.append(x_cg[j] / 27)
    x_cg_norm2.append(x_cg2[j] / 21)
s_h = []
s_h2 = []
l_t2 = []
l_t = []

for i in range(len(x_tail)):
    s_h.append(calculate_tail_area(x_cg_norm[i], x_tail[i], 27, 61, 2.2345, 0.05)[0])
    s_h2.append(calculate_tail_area(x_cg_norm2[i], x_tail[i], 21, 61, 2.2345, 0.05)[0])
    l_t.append(calculate_tail_area(x_cg_norm[i], x_tail[i], 27, 61, 2.2345, 0.05)[1])
    l_t2.append(calculate_tail_area(x_cg_norm2[i], x_tail[i], 21, 61, 2.2345, 0.05)[1])
plt.figure()
plt.plot(x_tail, [s_h * 0.80 for s_h in s_h], label='Fuselage length - 27 m', color='tab:blue')
plt.plot(x_tail, [s_h * 0.80 for s_h in s_h2], label='Fuselage length - 21 m', color='tab:green')

plt.xlabel(r'Normalized PE position ($x / l_{fuselage}$)')
plt.ylabel(f'Required horizontal tail area [$m^2$]')
plt.title('Required horizontal tail area with position of the PE')
plt.legend()
plt.ylim(0, 40)
plt.grid(True)

plt.figure()
plt.plot(x_tail, l_t, label='Fuselage length - 27 m', color='tab:blue')
plt.plot(x_tail, l_t2, label='Fuselage length - 21 m', color='tab:green')
plt.xlabel(r'Normalized PE position ($x / l_{fuselage}$)')
plt.ylabel(f'Tail arm length [$m$]')
plt.title('Tail arm length with PE position')
plt.legend()
plt.grid(True)

plt.show()


""" Test section"""

if __name__ == "__main__":
    length = np.linspace(0.2, 3, 150)
    cant_angles = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]

    plt.figure()

    # Color for cant = 30
    highlight_cant = 30
    highlight_color = "tab:blue"

    # Colormap for other cant lines
    cmap = plt.get_cmap("tab10")
    color_idx = 0

    for cant in cant_angles:
        weight = []
        for l in length:
            PE = PropulsiveEmpennage(l, cant)
            weight.append(PE.weight_ps())

        label = f"{cant}°"

        if cant == highlight_cant:
            plt.plot(length, weight, color=highlight_color, label=label)
        else:
            plt.plot(length, weight, color=cmap(color_idx % 10), linestyle='--', label=label)
            color_idx += 1

    plt.plot(1, 292.1, label="DUUC V0.1", marker='o', color='tab:blue')
    plt.xlabel(r'$Pylon\ length\ [m]$')
    plt.ylabel(r'$Mass\ [kg]$')
    plt.title(r'Pylon + Support Mass vs Pylon Length')
    plt.legend(title="Cant Angle")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

e_ratio = np.linspace(0.1, 2.0, 201)


def t_ratio(ratio):
    thrust = (2 * ratio) ** (1/3)
    return thrust

plt.figure()
plt.plot(e_ratio, [t_ratio(e) for e in e_ratio])
plt.axvline(x=1.1344, color="tab:blue", linestyle="dashed")
plt.axhline(y=t_ratio(1.1344), color="tab:blue", linestyle="dashed")
plt.plot(1.1344, t_ratio(1.1344), label="DUUC V0.1", marker='o', color='tab:blue')
plt.xlabel(r'Expansion ratio $\sigma_d$ [-]')
plt.ylabel(r'Thrust ratio between ducted and open rotor [-]')
plt.grid(True)
plt.title("Expansion ratio versus thrust ratio")
plt.legend()
plt.show()


from analysis_modules.factors import mass_of_section

class Duct:
    def __init__(self, geometry, geometry_pylon, geometry_support):
        super().__init__()
        self.duct_diameter = geometry[0]
        self.duct_chord = geometry[1]
        self.geometry_pylon = geometry_pylon
        self.geometry_support = geometry_support
        self.density = 0.54
    def weight(self):
        """ based on Torenbeek class II weight estimation for horizontal stabilizer """
        kh = 1.05
        sh = self.duct_diameter * np.pi * self.duct_chord
        vd = ref.v_dive
        sweep = np.radians(0)

        m_duct = kh * sh * (62 * (sh ** 0.2 * vd) / (1000 * np.sqrt(np.cos(sweep))) - 2.5)
        return m_duct

    def weight_ps(self):
        # Safety factors
        n_ult = 1.5  # Ultimate load factor
        k_stoot = 1.5  # Impact factor
        k_misc = 1.25  # Miscellaneous factor

        # Convert cant angle to radians
        c_rad = np.radians(self.geometry_pylon[3])

        # Total length of the beam
        l_tot = self.geometry_support[0] + self.geometry_pylon[0]  # [m]

        """ Material properties -> AL 6061-T6 """
        sigma_allow = 241 * 1e6  # Allowable stress [Pa] (converted from MPa to N/m^2)
        rho = 2700  # Mass per unit length of the beam [density]

        # Aerodynamic forces on pylon and support
        f_pylon = (0.5 * self.density * 128 ** 2 * 0.5 * self.geometry_pylon[1]
                   * self.geometry_pylon[0])  # [N]
        f_support = (0.5 * self.density * 128 ** 2 *0.5
                     * self.geometry_support[1] * self.geometry_support[0])  # [N]

        # Aerodynamic moment about the root
        m_aero = (0.5 * f_pylon * self.geometry_pylon[0] + (0.5 * self.geometry_support[0] + self.geometry_pylon[0]) * f_support) * np.cos(
            c_rad)  # [N·m]

        # Weight moment about the root
        m_weight = (self.weight() + 921 + 228 + 105) * (self.geometry_pylon[0] + self.geometry_support[0] * 0.5) * 9.81 * np.cos(
            c_rad)  # [N·m]

        # Sizing moment (maximum of aerodynamic and weight moments)
        m_sizing = max(m_aero, m_weight)  # [N·m]

        # Determine maximum thickness of the airfoil section
        num_list = [int(digit) for digit in self.geometry_pylon[2]]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness percentage
        thickness_pylon = (thickness / 100) * self.geometry_pylon[1]  # [m]

        num_list2 = [int(digit) for digit in self.geometry_support[2]]
        thickness_sup = num_list2[2] * 10 + num_list2[3]  # NACA thickness percentage
        thickness_support = (thickness_sup / 100) * self.geometry_support[1]  # [m]

        height = max(thickness_pylon, thickness_support) * 0.90  # Reduce maximum thickness by a safety factor [m]

        # Calculate width using bending stress formula
        width = (((6 * m_sizing * n_ult * k_stoot) / sigma_allow) / (height ** 2)) * k_misc  # [m]

        # Total mass of the beam
        mass_tot = height * width * l_tot * rho  # [kg]

        # Split mass into pylon and support sections
        m_pylon = mass_of_section(mass_tot, l_tot, 0, self.geometry_pylon[0])  # Mass of pylon section [kg]
        m_support = mass_of_section(mass_tot, l_tot, self.geometry_pylon[0], l_tot)  # Mass of support section [kg]

        return m_pylon, m_support


d_duct = np.linspace(0.5, 5, 151)

w_duct = []
w_support = []
sum_pe = []
for i in range(len(d_duct)):
    geometry = [d_duct[i], d_duct[i] / 2]
    geometry_pylon = [1, 1, "0012", 30]
    geometry_support = [d_duct[i], 0.75, "0012"]

    duct = Duct(geometry, geometry_pylon, geometry_support)
    w_duct.append(duct.weight())
    w_support.append(duct.weight_ps()[1])
    sum_pe.append(duct.weight() + duct.weight_ps()[1])


plt.figure()
plt.plot(d_duct, w_duct, label="Duct mass")
plt.plot(d_duct, w_support, label="Support mass")
#plt.plot(d_duct, [k + 921 + 228 + 105 for k in sum_pe], label="PE mass")
plt.xlabel(r"Duct diameter [m]")
plt.ylabel(r"Mass [kg]")
plt.title(r"Effect of duct diameter on PE mass")
plt.legend()
plt.grid(True)
plt.show()