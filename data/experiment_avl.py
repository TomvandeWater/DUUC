""" Ring wing and propulsive empennage modelled for power off, important to note that when scaling
the ring wing, the cl and cd need to be scaled with 10/AR to match the weissinger prediction model. """

a = [0, 3, 5, 7, 10, 12, 15]
cl_rw = [0, 1.388, 2.30911, 3.22262, 4.57318, 5.45647, 6.74936]
cd_rw = [0.00026, 0.01945, 0.05334, 0.10366, 0.20856, 0.29684, 0.45412]

cl_pe = [0.01453, 1.70731, 2.82823, 3.93870, 5.57683, 6.64542, 8.20478]
cd_pe = [0.00098, 0.03108, 0.08000, 0.15127, 0.29794, 0.42037, 0.63708]

cl_pe_nopylon = [-0.00408, 1.53547, 2.55691, 3.57030, 5.06774, 6.04606, 7.47564]
cd_pe_nopylon = [-0.00018, 0.02898, 0.07588, 0.14427, 0.28549, 0.40388, 0.61452]
