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
import Jeromy.analyzers.metallicity_separation   as metallicity_separation

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
from scipy.special import erfinv
from lmfit.models import GaussianModel
from tqdm import tqdm

def scale_spectrum(spectrum, exposure, bkg_reduction):
    yValues = spectrum.yValues
    xValues = spectrum.xValues

    yValues = [float(i) * float(exposure) / float(bkg_reduction) for i in yValues]
    xValues = [round(i, 2) for i in xValues]
    new_map = {}
    for i, x in enumerate(xValues):
        new_map[round(x, 2)] = yValues[i]
    spectrum.map = new_map
    spectrum.yValues = yValues
    spectrum.xValues = xValues

def gaussian(x, m, s, a):
    term1 = a / (s * np.sqrt(2*np.pi))
    term2 = (x-m)**2
    term3 = 2 * (s**2)
    return term1 * np.exp(-term2 / term3)

def orderOfMagnitude(number):
    return math.floor(math.log(number, 10))

if __name__ == "__main__":
    os.system("clear")
    args = create_parser()
    if args.configuration:
        config = configuration_loader.ConfigurationLoader(args.configuration)
    else:
        config = configuration_loader.ConfigurationLoader(systematics=args.systematics,
                                                          metallicity=args.metallicity,
                                                          exposure=args.exposure,
                                                          background_reduction=args.background_reduction)

        if args.metallicity == "high":
            systematics = args.systematics[:-2] + "lz"
            config2 = configuration_loader.ConfigurationLoader(systematics=systematics,
                                                               metallicity="low",
                                                               exposure=args.exposure,
                                                               background_reduction=args.background_reduction)
        else:
            systematics = args.systematics[:-2] + "hz"
            config2 = configuration_loader.ConfigurationLoader(systematics=systematics,
                                                               metallicity="high",
                                                               exposure=args.exposure,
                                                               background_reduction=args.background_reduction)
            
            
    syst = h.get_systematics(config)
    MetallicitySplit = metallicity_separation.MetallicitySplit(config, [1, 3], 1000)

    pred_cno_rate  = 250492.313050
    pred_cno_error = 124505.259759
    pred_cno_error = 116229.299903
    # pred_cno_error = 96762.552122

    lower_total  = MetallicitySplit.total_signal + pred_cno_rate - pred_cno_error
    middle_total = MetallicitySplit.total_signal + pred_cno_rate
    upper_total  = MetallicitySplit.total_signal + pred_cno_rate + pred_cno_error


    test_statistic_lower  = h.poisson_likelihood(MetallicitySplit.total_signal, lower_total )
    test_statistic_middle = h.poisson_likelihood(MetallicitySplit.total_signal, middle_total)
    test_statistic_upper  = h.poisson_likelihood(MetallicitySplit.total_signal, upper_total )

    h1_gaussian = h.gaussian(MetallicitySplit.test_likelihoods[0],
                             MetallicitySplit.test_fit.best_values["center"],
                             MetallicitySplit.test_fit.best_values["sigma"],
                             1)
    h0_gaussian = h.gaussian(MetallicitySplit.null_likelihoods[0],
                             MetallicitySplit.null_fit.best_values["center"],
                             MetallicitySplit.null_fit.best_values["sigma"],
                             1)

    max_value = max(max(h1_gaussian), max(h0_gaussian))

    
    plt.figure()
    plt.hist(MetallicitySplit.null_likelihoods[0], bins=MetallicitySplit.null_likelihoods[0],
             weights=MetallicitySplit.null_likelihoods[1], histtype="step")
    
    plt.plot(MetallicitySplit.null_likelihoods[0], h.gaussian(MetallicitySplit.null_likelihoods[0],
                                                              MetallicitySplit.null_fit.best_values["center"],
                                                              MetallicitySplit.null_fit.best_values["sigma"],
                                                              MetallicitySplit.null_fit.best_values["amplitude"]))
    
    plt.hist(MetallicitySplit.test_likelihoods[0], bins=MetallicitySplit.test_likelihoods[0],
             weights=MetallicitySplit.test_likelihoods[1], histtype="step")
    plt.plot(MetallicitySplit.test_likelihoods[0], h.gaussian(MetallicitySplit.test_likelihoods[0],
                                                              MetallicitySplit.test_fit.best_values["center"],
                                                              MetallicitySplit.test_fit.best_values["sigma"],
                                                              MetallicitySplit.test_fit.best_values["amplitude"]))

    plt.show()
    plt.figure()
    plt.plot(MetallicitySplit.null_likelihoods[0], h.gaussian(MetallicitySplit.null_likelihoods[0],
                                                              MetallicitySplit.null_fit.best_values["center"],
                                                              MetallicitySplit.null_fit.best_values["sigma"],
                                                              1))

    # plt.plot(MetallicitySplit.test_likelihoods[0], h.gaussian(MetallicitySplit.test_likelihoods[0],
    #                                                           MetallicitySplit.test_fit.best_values["center"],
    #                                                           MetallicitySplit.test_fit.best_values["sigma"],
    #                                                           1))
    plt.plot([test_statistic_lower , test_statistic_lower ], [0, max_value], "r--")
    plt.plot([test_statistic_middle, test_statistic_middle], [0, max_value], "r--")
    plt.plot([test_statistic_upper , test_statistic_upper ], [0, max_value], "r--")
    plt.ylim([0, max_value*1.2])
    plt.xlabel("Test statistic", fontsize=16)
    plt.ylabel("Likelihood", fontsize=16)
    plt.show()


    I, s = quad(h.gaussian, test_statistic_lower, np.inf, args=(MetallicitySplit.null_fit.best_values["center"],
                                                                MetallicitySplit.null_fit.best_values["sigma"],
                                                                1))
    CL = erfinv(1-I) * math.sqrt(2)
    
    
    print(f"p-value: {I}, CL: {CL}")
