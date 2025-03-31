""" Experimental investigation of lift, drag and pitching moment coefficient of five annular airfoils (October 1957)
    Source: https://ntrs.nasa.gov/citations/19930084906"""
#   Data extracted from graph reader for duct only
# cl vs alpha graphs (AR 1/3, 2/3, 1.0, 1.5, 3.0)

cla_a_ar_1_3 = [-4.129213483146068, -4.129213483146068, 0.6601123595505616, 0.6601123595505616, 4.522471910112358,
                8.539325842696627, 13.019662921348313, 17.345505617977526, 21.51685393258427, 21.51685393258427,
                25.379213483146067, 29.705056179775276, 33.87640449438202, 38.04775280898876, 42.37359550561797,
                52.57022471910112, 63.23033707865169, 73.42696629213482, 83.31460674157303, 94.28370786516854]
cla_cl_ar_1_3 = [-0.08014842300556588, -0.08014842300556588, -0.005936920222634509, -0.005936920222634509,
                 0.0712430426716141, 0.14545454545454548, 0.21669758812615958, 0.3057513914656772, 0.37699443413729133,
                 0.37699443413729133, 0.45417439703153994, 0.5461966604823748, 0.6352504638218924, 0.6174397031539889,
                 0.7183673469387756, 0.5788497217068647, 0.489795918367347, 0.3473098330241188, 0.14545454545454548,
                 0.06827458256029685]

cla_a_ar_2_3 = [-3.820224719101124, 0.19662921348314555, 4.213483146067414, 8.6938202247191, 13.019662921348313,
                17.191011235955052, 20.89887640449438, 29.396067415730336, 33.412921348314605, 37.89325842696629,
                42.52808988764045, 52.1067415730337, 62.61235955056179, 73.42696629213482, 83.77808988764045,
                94.28370786516854]
cla_cl_ar_2_3 = [-0.1424860853432282, 0, 0.13951762523191097, 0.27012987012987016, 0.4126159554730984, 0.549165120593692,
                 0.6886827458256031, 0.8549165120593692, 0.8460111317254175, 0.8816326530612246, 0.8846011131725419,
                 0.6797773654916512, 0.4482374768089054, 0.4690166975881262, 0.2077922077922078, 0.07421150278293136]

cla_a_ar_1_0 = [-3.974719101123596, 0.04213483146067354, 4.367977528089886, 8.6938202247191, 13.019662921348313,
                16.727528089887638, 21.05337078651685, 25.379213483146067, 28.932584269662918, 37.738764044943814,
                41.91011235955056, 52.57022471910112, 62.76685393258427, 73.5814606741573, 84.08707865168539,
                94.12921348314606]
cla_cl_ar_1_0 = [-0.2018552875695733, -0.0029684601113172545, 0.2018552875695733, 0.3977736549165121, 0.5936920222634509,
                 0.8252319109461967, 0.9766233766233767, 1.0419294990723562, 1.0003710575139149, 0.9439703153988869,
                 0.8935064935064936, 0.7213358070500928, 0.5432282003710576, 0.350278293135436, 0.2523191094619666,
                 0.08014842300556588]

cla_a_ar_1_5 = [-4.129213483146068, 0.04213483146067354, 4.213483146067414, 8.6938202247191, 12.71067415730337,
                17.036516853932582, 21.207865168539325, 25.07022471910112, 29.396067415730336, 33.412921348314605,
                37.58426966292134, 41.75561797752809, 52.57022471910112, 62.61235955056179, 73.5814606741573,
                83.4691011235955, 93.97471910112358]
cla_cl_ar_1_5 = [-0.3057513914656772, -0.011873840445269018, 0.2879406307977737, 0.5610389610389611, 0.8341372912801485,
                 1.1102040816326533, 1.3328385899814472, 1.1280148423005567, 1.0686456400742117, 1.0211502782931356,
                 0.9588126159554732, 0.9439703153988869, 0.8994434137291281, 0.6975881261595548, 0.4037105751391466,
                 0.18701298701298702, 0.08311688311688313]

cla_a_ar_3_0 = [-4.438202247191011, 0.04213483146067354, 4.67696629213483, 8.6938202247191, 12.556179775280896,
                17.191011235955052, 20.89887640449438, 25.07022471910112, 29.550561797752806, 37.738764044943814,
                42.06460674157303, 52.41573033707865, 62.92134831460673, 73.27247191011236, 83.93258426966291,
                93.82022471910112]
cla_cl_ar_3_0 = [-0.4304267161410019, -0.005936920222634509, 0.44230055658627093, 0.8430426716141003, 1.2259740259740262,
                 1.4129870129870132, 1.3476808905380335, 1.1933209647495362, 1.1844155844155846, 1.1695732838589983,
                 1.2289424860853433, 1.130983302411874, 0.9528756957328387, 0.6441558441558443, 0.29684601113172543,
                 0.0979591836734694]

# cd vs alpha graphs (AR 1/3, 2/3, 1.0, 1.5, 3.0)
cda_a_ar_1_3 = [-3.9290407358738504, 0.5519053876478317, 5.032851511169515, 9.369250985545335, 13.416557161629434,
                17.463863337713533, 21.944809461235216, 26.136662286465175, 30.039421813403415, 34.37582128777924,
                38.42312746386334, 42.47043363994744, 53.02233902759527, 63.429697766097235, 73.83705650459922,
                83.66622864651774, 94.65177398160316]
cda_cd_ar_1_3 = [0.01675977653631285, 0.01675977653631285, 0.025139664804469275, 0.03910614525139665,
                 0.06145251396648045, 0.0782122905027933, 0.111731843575419, 0.14245810055865923, 0.21508379888268156,
                 0.2905027932960894, 0.4329608938547486, 0.5698324022346369, 0.9888268156424581, 1.2039106145251397,
                 1.2150837988826817, 0.7793296089385475, 0.5614525139664804]

cda_a_ar_2_3 = [-4.073587385019711, 0.6964520367936924, 5.466491458607097, 9.658344283837057, 17.75295663600526,
                21.655716162943495, 25.992115637319316, 30.32851511169514, 34.37582128777924, 37.844940867279895,
                42.47043363994744, 52.73324572930355, 62.70696452036793, 73.69250985545335, 84.24441524310119,
                94.50722733245729]
cda_cd_ar_2_3 = [0.019553072625698324, 0.013966480446927375, 0.030726256983240226, 0.05027932960893855,
                 0.0893854748603352, 0.1452513966480447, 0.19273743016759778, 0.38268156424581007, 0.6312849162011174,
                 0.7486033519553073, 0.888268156424581, 1.1229050279329609, 1.1089385474860336, 1.2681564245810055,
                 0.9022346368715084, 0.547486033519553]

cda_a_ar_1_0 = [-3.639947437582129, 0.407358738501971, 5.032851511169515, 9.080157687253614, 13.561103810775297,
                17.897503285151117, 21.511169513797636, 25.992115637319316, 30.617608409986858, 33.942181340341655,
                38.42312746386334, 42.32588699080158, 53.31143232588699, 63.5742444152431, 73.98160315374507,
                83.8107752956636, 94.79632063074902]
cda_cd_ar_1_0 = [0.0223463687150838, 0.01675977653631285, 0.02793296089385475, 0.05027932960893855, 0.0782122905027933,
                 0.12849162011173185, 0.18156424581005587, 0.39385474860335196, 0.6340782122905028, 0.7653631284916201,
                 0.9189944134078213, 1.0251396648044693, 1.181564245810056, 1.3016759776536313, 1.223463687150838, 1,
                 0.5586592178770949]

cda_a_ar_1_5 = [-4.073587385019711, 0.11826544021024965, 9.080157687253614, 13.850197109067018, 17.75295663600526,
                22.23390275952694, 25.992115637319316, 30.473061760841, 34.37582128777924, 38.13403416557162,
                42.47043363994744, 53.16688567674113, 63.28515111695138, 73.98160315374507, 84.38896189224704,
                94.65177398160316]
cda_cd_ar_1_5 = [0.036312849162011177, 0.025139664804469275, 0.0558659217877095, 0.09217877094972067,
                 0.14245810055865923, 0.22625698324022347, 0.5139664804469274, 0.6312849162011174, 0.723463687150838,
                 0.8631284916201117, 0.9860335195530726, 1.1955307262569832, 1.4385474860335197, 1.3826815642458101,
                 0.8575418994413408, 0.5223463687150838]

cda_a_ar_3_0 = [-4.073587385019711, 0.6964520367936924, 4.888304862023654, 9.947437582128778, 13.416557161629434,
                18.04204993429698, 21.800262812089358, 26.136662286465175, 29.750328515111697, 38.71222076215506,
                43.04862023653088, 53.31143232588699, 63.28515111695138, 73.83705650459922, 84.09986859395532,
                94.79632063074902]
cda_cd_ar_3_0 = [0.04748603351955307, 0.03910614525139665, 0.0558659217877095, 0.08100558659217877, 0.10614525139664804,
                 0.2569832402234637, 0.38547486033519557, 0.6145251396648045, 0.7597765363128491, 1.0335195530726258,
                 1.1955307262569832, 1.5, 1.8324022346368716, 1.916201117318436, 1.2486033519553073, 0.7793296089385475]

# cm vs alpha graphs (AR 1/3, 2/3, 1.0, 1.5, 3.0)
cma_a_ar_1_3 = [-4.712643678160919, -0.22988505747126436, 4.022988505747127, 8.160919540229886, 12.183908045977011,
                16.666666666666668, 20.919540229885058, 25.057471264367816, 29.770114942528735, 33.793103448275865,
                38.04597701149425, 42.06896551724138, 52.758620689655174, 63.10344827586207, 73.79310344827586,
                83.79310344827586, 94.59770114942529]
cma_cm_ar_1_3 = [-0.02717391304347826, 0.002173913043478261, 0.035869565217391305, 0.057608695652173914,
                 0.08043478260869566, 0.09347826086956522, 0.1141304347826087, 0.12717391304347828, 0.14021739130434782,
                 0.14565217391304347, 0.10326086956521739, 0.06847826086956522, 0.01956521739130435,
                 -0.06847826086956522, -0.09565217391304348, -0.08260869565217391, -0.07934782608695652]

cma_a_ar_2_3 = [-4.597701149425287, -0.22988505747126436, 3.793103448275862, 8.275862068965518, 12.528735632183908,
                17.011494252873565, 21.149425287356323, 25.402298850574713, 29.54022988505747, 33.44827586206897,
                37.93103448275862, 42.52873563218391, 52.64367816091954, 63.2183908045977, 73.67816091954023,
                84.13793103448276, 94.71264367816092]
cma_cm_ar_2_3 = [-0.013043478260869566, 0.004347826086956522, 0.01847826086956522, 0.03695652173913044,
                 0.05652173913043478, 0.075, 0.0891304347826087, 0.10434782608695653, 0.021739130434782608,
                 -0.017391304347826087, -0.035869565217391305, -0.04673913043478261, -0.06304347826086956,
                 -0.07065217391304347, -0.1108695652173913, -0.11521739130434783, -0.08478260869565217]

cma_a_ar_1_0 = [-4.827586206896552, -0.3448275862068966, 3.793103448275862, 8.39080459770115, 12.528735632183908,
                17.011494252873565, 21.03448275862069, 25.517241379310345, 29.770114942528735, 33.9080459770115,
                38.39080459770115, 42.41379310344828, 52.87356321839081, 63.5632183908046, 73.44827586206897,
                84.02298850574712, 94.59770114942529]
cma_cm_ar_1_0 = [0.0032608695652173916, 0.011956521739130435, 0.02282608695652174, 0.029347826086956522,
                 0.041304347826086954, 0.059782608695652176, 0.08260869565217391, -0.07282608695652174,
                 -0.09891304347826087, -0.11304347826086956, -0.11521739130434783, -0.1141304347826087,
                 -0.10543478260869565, -0.09456521739130434, -0.08260869565217391, -0.09782608695652174,
                 -0.04782608695652174]

cma_a_ar_1_5 = [-4.827586206896552, -4.827586206896552, 0, 4.022988505747127, 8.045977011494253, 12.64367816091954,
                16.666666666666668, 20.919540229885058, 25.17241379310345, 29.42528735632184, 33.793103448275865,
                38.04597701149425, 42.298850574712645, 52.87356321839081, 63.333333333333336, 73.44827586206897,
                84.13793103448276, 94.82758620689656]
cma_cm_ar_1_5 = [0.008695652173913044, 0.008695652173913044, -0.0010869565217391304, 0, -0.002173913043478261,
                 0.005434782608695652, 0.025, 0.04782608695652174, -0.14891304347826087, -0.15869565217391304,
                 -0.15217391304347827, -0.14782608695652175, -0.1565217391304348, -0.14673913043478262, -0.15,
                 -0.12282608695652174, -0.08152173913043478, -0.0641304347826087]

cma_a_ar_3_0 = [-4.482758620689655, -0.22988505747126436, 4.252873563218391, 8.620689655172413, 12.873563218390805,
                16.551724137931036, 21.60919540229885, 25.517241379310345, 30, 34.02298850574713, 38.160919540229884,
                42.758620689655174, 53.2183908045977, 63.79310344827586, 73.9080459770115, 84.36781609195403,
                94.71264367816092]
cma_cm_ar_3_0 = [0.035869565217391305, -0.010869565217391304, -0.058695652173913045, -0.09782608695652174,
                 -0.10434782608695653, -0.29130434782608694, -0.3391304347826087, -0.21630434782608696,
                 -0.19782608695652174, -0.20217391304347826, -0.21847826086956523, -0.2423913043478261,
                 -0.2597826086956522, -0.2902173913043478, -0.2619565217391304, -0.15108695652173912,
                 -0.08586956521739131]

# cl vs cd graphs (AR 1/3, 2/3, 1.0, 1.5, 3.0)
clcd_cd_1_3 = [0.026192170818505337, 0.021637010676156584, 0.02277580071174377, 0.029608540925266904,
               0.05921708185053381, 0.0899644128113879, 0.10021352313167259, 0.1412099644128114,
               0.20384341637010675, 0.2858362989323843]
clcd_cl_1_3 = [-0.07428571428571429, 0, 0.06857142857142857, 0.14857142857142858, 0.22285714285714286,
               0.30857142857142855, 0.37714285714285717, 0.45714285714285713, 0.5485714285714286, 0.64]

clcd_cd_2_3 = [0.024858757062146894, 0.015819209039548022, 0.021468926553672316, 0.04067796610169492,
               0.07344632768361582, 0.0847457627118644, 0.13333333333333333, 0.19322033898305085]
clcd_cl_2_3 = [-0.14285714285714285, -0.011428571428571429, 0.13714285714285715, 0.2742857142857143,
               0.4114285714285714, 0.5485714285714286, 0.6914285714285714, 0.8171428571428572]

clcd_cd_1_0 = [0.02272727272727273, 0.010227272727272729, 0.01931818181818182, 0.05795454545454546,
               0.08863636363636365, 0.12272727272727274, 0.17272727272727276]
clcd_cl_1_0 = [-0.2, 0, 0.2057142857142857, 0.4114285714285714, 0.6057142857142858, 0.84, 0.9885714285714285]

clcd_cd_1_5 = [0.028436018957345967, 0.026161137440758292, 0.026161137440758292, 0.052322274881516584,
               0.07393364928909951, 0.13649289099526066, 0.22407582938388623]
clcd_cl_1_5 = [-0.28, 0.017142857142857144, 0.30857142857142855, 0.5771428571428572, 0.8571428571428571,
               1.12, 1.3428571428571427]

clcd_cd_3_0 = [0.06978723404255319, 0.04765957446808511, 0.06638297872340425, 0.1072340425531915, 0.14468085106382977]
clcd_cl_3_0 = [-0.4228571428571429, 0.017142857142857144, 0.46285714285714286, 0.8742857142857143, 1.262857142857143]
