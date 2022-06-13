#!/usr/bin/env python3

import os
import math
import numpy as np
import matplotlib.pyplot as plt
import Jeromy
import Jeromy.statistics.test_statistic as test_statistic
import Jeromy.statistics.test_ensemble  as test_ensemble
import Jeromy.statistics.run_chi2_tests as run_chi2_tests

import Jeromy.core.true_solar               as true_solar
import Jeromy.analyzers.counting_experiment as counting_experiment
import Jeromy.analyzers.likelihood_ratios   as likelihood_ratios

import Jeromy.IO.FileImporter               as FI
import Jeromy.IO.ImportDataHighZ            as DataHZ
import Jeromy.IO.ImportDataLowZ             as DataLZ
import Jeromy.ToyMCGenerator.ToyMonteCarlo  as ToyMonteCarlo
import Jeromy.IO.UncertaintyImporter        as UncIm
import Jeromy.config.configuration_loader   as configuration_loader
import Jeromy.plots.make_spectrum_plots     as make_spectrum_plots
import Jeromy.plots.make_true_plots         as make_true_plots
import Jeromy.utilities.helpers             as h

from Jeromy.utilities.argument_parser import create_parser
from scipy.integrate import quad
from lmfit.models import GaussianModel
from tqdm import tqdm
import yaml

def gaussian(x, m, s, a):
    term1 = a / (s * np.sqrt(2*np.pi))
    term2 = (x-m)**2
    term3 = 2 * (s**2)
    return term1 * np.exp(-term2 / term3)

with open(Jeromy.__jeromy_data__+"/solar_rates.yml", "r") as stream:
    solar_rates = yaml.safe_load(stream)


if __name__ == "__main__":
    os.system("clear")
    args = create_parser()
    if args.configuration:
        config = configuration_loader.ConfigurationLoader(args.configuration)
        if config.metallicity == "high": temp_metallicity = "low"
        if config.metallicity == "low": temp_metallicity = "high"
        
        config2 = configuration_loader.ConfigurationLoader(systematics=config.systematic_uncertainties.name,
                                                           metallicity=temp_metallicity,
                                                           exposure=config.exposure,
                                                           background_reduction=config.background_reduction)
    else:
        config = configuration_loader.ConfigurationLoader(systematics=args.systematics,
                                                          metallicity=args.metallicity,
                                                          exposure=args.exposure,
                                                          background_reduction=args.background_reduction)
        if args.metallicity == "high": temp_metallicity = "low"
        if args.metallicity == "low": temp_metallicity = "high"
        
        config2 = configuration_loader.ConfigurationLoader(systematics=args.systematics,
                                                           metallicity=temp_metallicity,
                                                           exposure=args.exposure,
                                                           background_reduction=args.background_reduction)

    systematics = config.systematic_uncertainties
    


    O15_rate = solar_rates[config.metallicity]["O15"]/10
    F17_rate = solar_rates[config.metallicity]["F17"]/10
    true_solar_rate = (O15_rate +  F17_rate)
    O15_error = [systematics["O15_flux_e"], systematics["Sur_prob_e"], systematics["trigger_e"]]
    F17_error = [systematics["F17_flux_e"], systematics["Sur_prob_e"], systematics["trigger_e"]]
    O15_error = h.add_in_quadrature(O15_error)*O15_rate
    F17_error = h.add_in_quadrature(F17_error)*F17_rate
    CNO_error = h.add_in_quadrature([O15_error, F17_error])

    test = counting_experiment.CountingExperiment("H0", [1, 3], config, True)

    x_values = np.linspace(true_solar_rate-5*CNO_error, true_solar_rate+5*CNO_error, 50)
    y_values = [gaussian(i, true_solar_rate, CNO_error, 1) for i in x_values]
    c_values = [gaussian(i, test.cno_rate, test.cno_uncertainty, 1) for i in x_values]

    print(f"Prediction:  {true_solar_rate} ± {CNO_error}")
    print(f"Calculation: {round(test.cno_rate, 2)} ± {round(test.cno_uncertainty, 2)}")

    q = 0
    for i in range(len(y_values)):
        if y_values[i]:
            q += (c_values[i] - y_values[i])**2 / y_values[i]

    print(q)
    plt.figure()
    plt.plot(x_values, y_values, label="prediction")
    plt.plot(x_values, c_values, label="calculation")
    plt.legend(loc="best")
    plt.show()
