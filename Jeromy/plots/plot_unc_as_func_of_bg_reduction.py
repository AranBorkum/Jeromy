import numpy as np
import matplotlib.pyplot as plt

directory = "/Users/aranborkum/Docterate/Projects/SolarStuff/CNOFluxAnalysis/Jeromy/"

standard1000 = np.loadtxt(directory+"StandardTo1000.txt"         , delimiter=",", dtype=str)
standard     = np.loadtxt(directory+"CNOdiscoveryStandard.txt"   , delimiter=",", dtype=str)
optimistic   = np.loadtxt(directory+"CNOdiscoveryOptimistic.txt" , delimiter=",", dtype=str)
betterBoron  = np.loadtxt(directory+"CNOdiscoveryBetterBoron.txt", delimiter=",", dtype=str)
bestBoron    = np.loadtxt(directory+"CNOdiscoveryBestBoron.txt"  , delimiter=",", dtype=str)
perfect      = np.loadtxt(directory+"CNOdiscoveryPerfect.txt"    , delimiter=",", dtype=str)

standard_x = [int(i) for i in standard[:,  1]]
standard_y = []
for i in range(len(standard[:, 1])):
    standard_y.append(100* float(standard[:, -2][i]) / float(standard[:, -3][i]))

optimistic_x = [int(i) for i in optimistic[:,  1]]
optimistic_y = []
for i in range(len(optimistic[:, 1])):
    optimistic_y.append(100* float(optimistic[:, -2][i]) / float(optimistic[:, -3][i]))

betterBoron_x = [int(i) for i in betterBoron[:,  1]]
betterBoron_y = []
for i in range(len(betterBoron[:, 1])):
    betterBoron_y.append(100* float(betterBoron[:, -2][i]) / float(betterBoron[:, -3][i]))

bestBoron_x = [int(i) for i in bestBoron[:,  1]]
bestBoron_y = []
for i in range(len(bestBoron[:, 1])):
    bestBoron_y.append(100* float(bestBoron[:, -2][i]) / float(bestBoron[:, -3][i]))

perfect_x = [int(i) for i in perfect[:,  1]]
perfect_y = []
for i in range(len(perfect[:, 1])):
    perfect_y.append(100* float(perfect[:, -2][i]) / float(perfect[:, -3][i]))

standard1000_x = [int(i) for i in standard1000[:,  1]]
standard1000_y = []
for i in range(len(standard1000[:, 1])):
    standard1000_y.append(100* float(standard1000[:, -2][i]) / float(standard1000[:, -3][i]))

    
plt.figure()
# plt.plot(standard_x   , standard_y   , label="standard")
plt.plot(standard1000_x   , standard1000_y   , label="standard", linewidth=2)
# plt.plot(optimistic_x , optimistic_y , label="optimistic")
# plt.plot(betterBoron_x, betterBoron_y, label="better boron")
# plt.plot(bestBoron_x  , bestBoron_y  , label="best boron")
# plt.plot(perfect_x    , perfect_y    , label="perfect")
plt.ylabel("Uncertainty as percentage of calculated CNO rate", fontsize=24)
plt.xlabel("Background reduction factor", fontsize=24)
# plt.legend(loc="best")
plt.yscale("log")
plt.xscale("log")
plt.tick_params(axis='both', which='major', labelsize=20)
plt.show()
