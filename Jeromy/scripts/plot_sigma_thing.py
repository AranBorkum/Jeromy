import numpy as np
import matplotlib.pyplot as plt

lowZ = np.loadtxt("/Users/aranborkum/Docterate/Projects/SolarStuff/CNOFluxAnalysis/Jeromy/standard_unc_lowZ_varBkgRed.txt", delimiter=",", dtype=str)
highZ = np.loadtxt("/Users/aranborkum/Docterate/Projects/SolarStuff/CNOFluxAnalysis/Jeromy/standard_unc_highZ_varBkgRed.txt", delimiter=",", dtype=str)

plt.figure()
plt.plot(lowZ[:, 1], lowZ[:, 4])
plt.plot(highZ[:, 1], highZ[:, 4])
plt.xscale("log")
plt.show()
