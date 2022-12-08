from energy_loss_3 import *
from plots import plot_histogram
from histogram_fitters import *
print ("beginning testing")
data = energy_data()
chx = [0.00000000e+00, 6.66666667e-05, 1.33333333e-04, 2.00000000e-04,
 2.66666667e-04, 3.33333333e-04, 4.00000000e-04, 4.66666667e-04,
 5.33333333e-04, 6.00000000e-04]
chy = [1500.6320261106, 1426.7784637038178, 1351.7497554118797, 1276.2285855547163, 1189.0, 1100.2522614231375, 1009.0, 912.9896904539689, 813.6829992511642, 711.4027979465324]   
chmu = [1467.4529594102585, 1393.2297107781342, 1314.273037553567, 1230.6085276000981, 1142.2632945927603, 1049.2663979892923, 951.6490994383291, 849.4472616946429, 742.7245010393489, 631.7073718233077]


chx = [0.0, 1.1999999999999999e-05, 2.3999999999999997e-05, 3.5999999999999994e-05, 4.7999999999999994e-05, 5.9999999999999995e-05, 7.199999999999999e-05, 8.4e-05, 9.599999999999999e-05, 0.00010799999999999998, 0.00011999999999999999, 0.00013199999999999998, 0.000168, 0.00020399999999999997, 0.00021599999999999996, 0.00022799999999999999, 0.00025199999999999995, 0.00026399999999999997, 0.000276, 0.00028799999999999995, 0.0003, 0.00031199999999999994, 0.00032399999999999996, 0.000336, 0.00034799999999999995, 0.00035999999999999997, 0.00037199999999999993, 0.00038399999999999996, 0.000396, 0.00040799999999999994, 0.00041999999999999996, 0.00043199999999999993, 0.00044399999999999995, 0.00045599999999999997, 0.00046799999999999994, 0.00047999999999999996, 0.0004919999999999999, 0.0005039999999999999, 0.000516, 0.0005279999999999999, 0.0005399999999999999, 0.000552, 0.0005639999999999999, 0.0005759999999999999, 0.000588, 0.0006]
chy = [1507.918868094196, 1494.4041425025653, 1480.7632131418038, 1466.9950050723264, 1453.1000767158066, 1439.0797240835786, 1424.9361604406115, 1410.6725609218033, 1396.2904144609938, 1385.6594209250368, 1370.429709693116, 1356.4226609484572, 1315.0, 1269.0, 1248.850512716216, 1239.0, 1207.0, 1189.839433321339, 1176.0, 1160.0, 1144.0, 1127.4628815464494, 1111.1347564152065, 1095.6163418968154, 1074.7034145694686, 1062.0882686952016, 1046.6351117445047, 1030.0, 1011.9331209431034, 995.8659288957315, 979.0, 960.663473661209, 944.7080223265712, 928.0, 904.7842569564281, 893.4804783739763, 875.8178717363379, 856.91908894877, 840.0256142581468, 821.3958293325709, 804.9034660161163, 785.1037006081465, 769.3952244115728, 751.0184013642908, 732.0, 712.5906000104804]
chmu = [1474.7510425514602, 1461.1929102062534, 1447.4972574832454, 1433.6646845615917, 1419.6954384886371, 1405.5896542350297, 1391.3474611822185, 1376.9689903087274, 1362.4543742285334, 1347.8037471938405, 1333.0172450982845, 1318.0950054806551, 1272.5152616569076, 1225.7189905398382, 1209.850578678853, 1193.8475895530792, 1161.4384919456643, 1145.0326940877908, 1128.492940273476, 1111.8193903713368, 1095.012206164442, 1078.071551430875, 1060.9975920588267, 1043.7904962173438, 1026.450434618268, 1008.9775809282205, 991.3721124259172, 973.6342110549487, 955.7640651011598, 937.7618718328022, 919.6278415850484, 901.362203949952, 882.9652169456689, 864.4371802759001, 845.7784540364127, 826.9894844575433, 808.0708384615114, 789.0232489294942, 769.8476725842452, 750.5453622721024, 731.1179551542427, 711.5675778832642, 691.896969253627, 672.1096200936197, 652.2099293461192, 632.2033744097957]

chx = [0.0, 1.1999999999999999e-05, 2.3999999999999997e-05, 3.5999999999999994e-05, 4.7999999999999994e-05, 5.9999999999999995e-05, 7.199999999999999e-05, 8.4e-05, 9.599999999999999e-05, 0.00010799999999999998, 0.00011999999999999999, 0.00013199999999999998, 0.00014399999999999998, 0.00015599999999999997, 0.000168, 0.00017999999999999998, 0.00019199999999999998, 0.00020399999999999997, 0.00021599999999999996, 0.00022799999999999999, 0.00025199999999999995, 0.00026399999999999997, 0.000276, 0.00028799999999999995, 0.0003, 0.00031199999999999994, 0.00032399999999999996, 0.000336, 0.00034799999999999995, 0.00035999999999999997, 0.00037199999999999993, 0.00038399999999999996, 0.000396, 0.00040799999999999994, 0.00041999999999999996, 0.00043199999999999993, 0.00044399999999999995, 0.00045599999999999997, 0.00046799999999999994, 0.00047999999999999996, 0.0004919999999999999, 0.0005039999999999999, 0.000516, 0.0005279999999999999, 0.0005399999999999999, 0.000552, 0.0005639999999999999, 0.000588, 0.0006]
chy = [1523.6722879892823, 1508.7778217601676, 1493.7803560101654, 1478.6874098123103, 1463.6717605826866, 1448.254819083511, 1432.936647171258, 1417.7258480043852, 1402.6700450562084, 1387.7508455021655, 1372.9813629936104, 1356.9045075686051, 1342.4485305837911, 1326.6857091099494, 1312.5817962476224, 1295.5694534113811, 1286.7156841177352, 1269.0, 1253.213345213271, 1237.6853749751108, 1207.0, 1191.3733959086828, 1175.0, 1159.9502780576179, 1142.8054971347112, 1127.0, 1111.0, 1095.5484457188518, 1079.402179606987, 1062.0, 1046.211290013742, 1027.5861935750893, 1011.225213250472, 995.099879491232, 979.9329984979094, 962.0, 945.0, 926.8432824405935, 903.033394827316, 891.8086634908855, 871.7535563321991, 858.0, 840.0, 823.0, 806.2536469230793, 787.0, 769.0, 732.260902103773, 712.4211270728831]
chmu = [1490.9660074063431, 1476.1040211271918, 1461.1481752054317, 1446.09865895771, 1430.9555670970294, 1415.7189949720162, 1400.3890385685727, 1384.9657945121999, 1369.44936007067, 1353.8398331571118, 1338.1373123335907, 1322.3418968152698, 1306.4536864752695, 1290.4727818503632, 1274.3992841476866, 1258.2332952526806, 1241.9749177385486, 1225.6242548775986, 1209.1814106549605, 1192.6464897853598, 1159.3008407426894, 1142.49032586477, 1125.5881610109652, 1108.5944550145637, 1091.5093177254835, 1074.33286015047, 1057.065194665242, 1039.7064353368855, 1022.2566984111402, 1004.7161030399122, 987.0847723494131, 969.3628349786302, 951.5504272508487, 933.6476961769185, 915.6548035268916, 897.571931245373, 879.3992885242526, 861.1371208831418, 842.7857216417506, 824.3454461985723, 805.8167295558302, 787.2001075510464, 768.4962422704837, 749.705952128878, 730.8302471033543, 711.8703696073543, 692.8278414830772, 654.5026463592947, 635.2249379806854]

#energies = get_energies(data, verbose = True, plot = True)
#thickness = get_thickness(energies, verbose = True)
#do_gold_thickness(thickness)
#a = 5 / 0
#thiccness_dx('gold', energies['empty'].val, energies[('gold')].val)
#thiccness_dx('2gold', energies['empty'].val, energies[('2gold')].val)
#thiccness_dx('4gold', energies['empty'].val, energies[('4gold')].val, 0.34490747, 0.16126968)
metadata = ('4gold', 0, 1)
#plot_histogram(metadata, data[metadata][1])
#histogram = data[metadata][1]
#result = fit_histogram(histogram, const = True, soft = True, plot = False)
#exiting = optimal_energy_fit(data[metadata][1])
energies, e = get_energies(data, channel_map=True)
ey = [e(ch).val for ch in chy]
emu = [e(ch).val for ch in chmu]
incident_mu = np.dot(data[('empty', 0, 1)][1], np.array(range(2048)) / np.sum(data[('empty', 0, 1)][1]))
i_e_mu = e(incident_mu).val
print (i_e_mu)
def thiccness_dx_2(foil, incident, incident2, exiting, c_a = c_a['gold'].val, c_b = c_b['gold'].val):
    print (incident)
    E = incident + m_a
    x = 0
    dt = 2e-16
    Es = [incident]
    xs = [x]
    #dEs = [dEdx_simple(b(E), foil) * v(b(E)) * 100 * dt]
    
    Cs = []
    L1s = []
    L2s = []
    Las = []
    i = 0
    E = incident + m_a
    dEs = [dEdx(b(E), foil, c_a, c_b)]
    x = 0
    dx = 1e-16 * v(b(E)) * 100
    #steps = 0
    while E > exiting + m_a:
        beta = b(E)
        x += dx
        dE = dEdx(beta, foil, c_a, c_b) * dx
        E-= dE
        Es.append(E - m_a)
        i += 1
        dEs.append(dE / dx)
        Las.append(K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * density_map[foil] * L_a(beta, foil))
        Cs.append(K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * density_map[foil] * -C_Z(beta, foil, c_a, c_b))
        L1s.append(K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * density_map[foil] * 2*L1(beta, foil, c_a, c_b))
        L2s.append(K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * density_map[foil] * L2(beta))
        xs.append(x)
        #if i % 50 == 0:
        #    plt.plot(xs, Es, label = "Es")
        #    plt.plot(xs, [d*20 for d in dEs], label = "dEs")
        #    plt.axhline(exiting, color = 'r', ls = '--')
        #    plt.legend()
        #    plt.show()
    #print (x)

    print (incident2)
    E2 = incident2 + m_a
    x2 = 0
    dt2 = 2e-16
    Es2 = [incident2]
    xs2 = [x2]
    #dEs = [dEdx_simple(b(E), foil) * v(b(E)) * 100 * dt]
    
    Cs2 = []
    L1s2 = []
    L2s2 = []
    Las2 = []
    i = 0
    E2 = incident2 + m_a
    dEs2 = [dEdx(b(E), foil, c_a, c_b)]
    x2 = 0
    dx2 = 1e-16 * v(b(E)) * 100
    #steps = 0
    while E2 > exiting + m_a:
        beta = b(E2)
        x2 += dx2
        dE2 = dEdx(beta, foil, c_a, c_b) * dx2
        E2-= dE2
        Es2.append(E2 - m_a)
        i += 1
        dEs2.append(dE2 / dx2)
        Las2.append(K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * density_map[foil] * L_a(beta, foil))
        Cs2.append(K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * density_map[foil] * -C_Z(beta, foil, c_a, c_b))
        L1s2.append(K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * density_map[foil] * 2*L1(beta, foil, c_a, c_b))
        L2s2.append(K * 2**2 * Z_map[foil] / mass_map[foil] / beta**2 * density_map[foil] * L2(beta))
        xs2.append(x2)





    plt.title("E(x) integration")
    plt.plot(xs, Es, label = "E(x) starting from ⟨E_α⟩", color = 'purple')
    plt.plot(xs2, Es2, label = "E(x) starting from ⟨E_α(0)⟩ calculated from histogram", color = 'magenta')
    plt.axhline(exiting, color = 'r', ls = '--', label = "3x gold E_exiting")
    plt.axhline(energies['2gold'].val, color = 'g', ls = '--', label = "2x gold E_exiting")
    plt.axhline(energies['gold'].val, color = 'b', ls = '--', label = "1x gold E_exiting")
    #plt.plot(chx, ey, ls = 'none', marker = '.', label = "Moyal center prediction", color = 'orange')
    plt.plot(chx, emu, ls = 'none', marker = '.', label = "Moyal expectation prediction", color = 'black')
    plt.xlabel("x (cm)")
    plt.ylabel("Kinetic Energy (MeV)")
    plt.ylim([0, 6])
    plt.legend()
    plt.savefig("E_integration_with_prediction.png")
    plt.show()
    #
    #plt.plot(xs, [d*20 for d in dEs], label = "dEs*20")
    #plt.plot(xs, dE_complexs, label = "With 'small' corrections")
    plt.title("dE/dx integration")
    plt.plot(xs, dEs, label = "-dE/dx", color = "red")
    labels = ["C/Z contribution", "L_1  contribution", "L_2  contribution", "L_a contribution"]
    cols = ["magenta", "orange", "green", "blue"]
    for i,y in enumerate((Cs, L1s, L2s, Las)):
        plt.plot(xs[1:], y, label = labels[i], color = cols[i])
    plt.xlabel("x (cm)")
    plt.ylabel("Kinetic Energy Derivative (MeV/cm)")
    plt.legend()
    plt.savefig("dE_dx_integration.png")
    plt.show()
    #print (xs)
    #print (x)
    #print (len(xs))
    #print (dEs)
    #return x, interpolate(Es, xs)[0]
    #return xs, Es  
"""
metadata = ('2gold', 0, 1)
plot_histogram(metadata, data[metadata][1])
print (optimal_energy_fit(data[metadata][1], plot = False))
metadata = ('4gold', 0, 1)
plot_histogram(metadata, data[metadata][1])
print (optimal_energy_fit(data[metadata][1], plot = False))
metadata = ('gold', 0, 1)
plot_histogram(metadata, data[metadata][1])
print (optimal_energy_fit(data[metadata][1], plot = False))
metadata = ('2gold', 0, 1)
plot_histogram(metadata, data[metadata][1])
print (optimal_energy_fit(data[metadata][1], plot = False))
metadata = ('4gold', 0, 1)
plot_histogram(metadata, data[metadata][1])
print (optimal_energy_fit(data[metadata][1], plot = False)
"""
#redchi = fit_histogram(data[('empty',0,1)][1], min_counts = 13, plot = True, foil = 'empty')[1]

#for j in range(5, 16):
#    redchi = fit_histogram(data[('empty',0,1)][1], min_counts = j, plot = False, foil = 'empty')[1]
#    print ("j:",j,"| Redchi:",redchi)
#print (exiting.val)
thiccness_dx_2(metadata[0], e_empty.val, i_e_mu, energies[metadata[0]].val)
#for incident in (e_empty.val, i_e_mu):
#    thiccness_dx_2(metadata[0], incident, energies[metadata[0]].val)
