
import yaml
import Jeromy
import Jeromy.utilities.helpers as h
import numpy as np
import matplotlib.pyplot as plt

with open(Jeromy.__jeromy_data__+"/solar_rates.yml", "r") as stream:
    solar_rates = yaml.safe_load(stream)

def gaussian(x, m, s, a):
    term1 = a / (s * np.sqrt(2*np.pi))
    term2 = (x-m)**2
    term3 = 2 * (s**2)
    return term1 * np.exp(-term2 / term3)
    
def make_true_rate_plots(config):

    systematics = config.systematic_uncertainties
    print(solar_rates)
    O15_rate = solar_rates[config.metallicity]["O15"]
    F17_rate = solar_rates[config.metallicity]["F17"]
    CNO_rate = O15_rate + F17_rate
    
    O15_error = [systematics["O15_flux_e"], systematics["Sur_prob_e"], systematics["trigger_e"]]
    F17_error = [systematics["F17_flux_e"], systematics["Sur_prob_e"], systematics["trigger_e"]]
    O15_error = h.add_in_quadrature(O15_error)*O15_rate
    F17_error = h.add_in_quadrature(F17_error)*F17_rate
    CNO_error = h.add_in_quadrature([O15_error, F17_error])

    O15_xValues = np.linspace(O15_rate-(O15_error*7), O15_rate+(O15_error*7), 1000)
    O15_yValues = [gaussian(i, O15_rate, O15_error, 1) for i in O15_xValues]
       
    F17_xValues = np.linspace(F17_rate-(F17_error*7), F17_rate+(F17_error*7), 1000)
    F17_yValues = [gaussian(i, F17_rate, F17_error, 1) for i in F17_xValues]

    CNO_xValues = np.linspace(CNO_rate-(CNO_error*7), CNO_rate+(CNO_error*7), 1000)
    CNO_yValues = [gaussian(i, CNO_rate, CNO_error, 1) for i in CNO_xValues]

    plt.figure()
    plt.subplot(1, 3, 1)
    label = r"Mean: %.2E, $\sigma$: %.2E (%i%%)"%(O15_rate, O15_error, 100*O15_error/O15_rate)
    plt.plot(O15_xValues, O15_yValues, label=label)
    plt.legend(loc="best")
    
    plt.subplot(1, 3, 2)
    label = r"Mean: %.2E, $\sigma$: %.2E (%i%%)"%(F17_rate, F17_error, 100*F17_error/F17_rate)
    plt.plot(F17_xValues, F17_yValues, label=label)
    plt.legend(loc="best")
    
    plt.subplot(1, 3, 3)
    label = r"Mean: %.2E, $\sigma$: %.2E (%i%%)"%(CNO_rate, CNO_error, 100*CNO_error/CNO_rate)
    plt.plot(CNO_xValues, CNO_yValues, label=label)
    plt.legend(loc="best")
    plt.show()
