#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import Jeromy
import Jeromy.statistics.test_statistic as test_statistic
import Jeromy.statistics.test_ensemble  as test_ensemble
import Jeromy.statistics.run_chi2_tests as run_chi2_tests

import Jeromy.core.true_solar               as true_solar
import Jeromy.analyzers.counting_experiment as counting_experient
import Jeromy.analyzers.likelihood_ratios   as likelihood_ratios

import Jeromy.ToyMCGenerator.ToyMonteCarlo  as ToyMonteCarlo
import Jeromy.IO.UncertaintyImporter        as UncIm
import Jeromy.config.configuration_loader   as configuration_loader
import Jeromy.plots.make_spectrum_plots     as make_spectrum_plots
import Jeromy.plots.make_true_plots         as make_true_plots
import Jeromy.utilities.helpers             as h

from Jeromy.utilities.argument_parser import create_parser
from scipy.integrate import quad
from lmfit.models import GaussianModel


def gaussian(x, m, s, a):
    term1 = a / (s * np.sqrt(2*np.pi))
    term2 = (x-m)**2
    term3 = 2 * (s**2)
    return term1 * np.exp(-term2 / term3)

if __name__ == "__main__":

    os.system("clear")
    args = create_parser()

    if args.config_help: configuration_loader.dump_configurations()

    else:
        config = configuration_loader.ConfigurationLoader(args.configuration)
        config.print()
        if args.exposure: config.set_exposure(args.exposure)
        if args.background_reduction: config.set_background_reduction(args.background_reduction)
        if args.systematics: config.set_systematic_uncertainties(args.systematics)
        if args.metallicity: config.set_metallicity(args.metallicity)


        # A helper function to dump all of the systematics used within the analysis.
        if args.dump_systematics: h.dump_systematics(config)

        # Plot dumper for various outputs. These include the various spectra that
        # go into the analysis and the probabity distributions for the CNO neutinos,
        # individually and combined.
        if args.make_plots:
            make_spectrum_plots.make_spectra(config)
            # make_true_plots.make_true_rate_plots(config)

        # test_ensemble.test_ensemble(config)
        test = counting_experient.CountingExperiment("Daniel", [1, 4], config, False, True)
        print(f"Configuration: {config.name}\n"
              "----------------------------------------------------\n"
              "Rate: %f Uncertainty: %f (%i%%)\n" %
              (test.cno_rate, test.cno_uncertainty, test.cno_uncertainty/test.cno_rate*100))
        testb8 = counting_experient.CountingExperiment_B8("Daniel", [5, 8], config, False)
        print(f"Configuration: {config.name}\n"
              "----------------------------------------------------\n"
              "Rate: %.2E Uncertainty: %.2E (%i%%)\n" %
              (testb8.b8_rate, testb8.b8_uncertainty, testb8.b8_uncertainty/testb8.b8_rate*100))

"""
        
        likelihood = likelihood_ratios.LikelihoodTest(test, test_statistic.BakerCousinsChi, config, False)
        print(likelihood._predicted.mean, likelihood._predicted.error)
        path = f"{Jeromy.__jeromy_output__}/"
        path += f"{config.metallicity}Z_"
        path += f"exp{float(config.exposure) * 10}kt_year_"
        path += f"{args.systematics}"
        h.ensure_directory_exists(path)
        plt.figure()
        plt.subplot(1, 3, 1)
        plt.hist(likelihood.bins, bins=likelihood.bins, weights=likelihood.observed , label="obs", histtype="step")
        plt.title("Observed")
        plt.subplot(1, 3, 2)
        plt.hist(likelihood.true_bins, bins=likelihood.true_bins, weights=likelihood.true_pred, label="pre",
                 histtype="step")
        plt.title("Predicted")
        plt.subplot(1, 3, 3)
        plt.hist(likelihood.bins, bins=likelihood.bins, weights=likelihood.observed , label="obs", histtype="step")
        plt.hist(likelihood.bins, bins=likelihood.bins, weights=likelihood.predicted, label="pre", histtype="step")
        plt.legend()
        plt.savefig(f"{path}/likelihoods.png")
        if args.make_plots: plt.show()
        plt.close()
        
        print(likelihood.chi2)
        likelihood_ensemble1 = likelihood_ratios.LikelihoodTestEnsemble("Jim", [1, 3], config, 1000, False)
        likelihood_ensemble2 = likelihood_ratios.LikelihoodTestEnsemble("Jim", [1, 3], config, 1000, True )

        plt.figure()
        plt.subplot(1, 2, 1)
        hist1, fit1 = likelihood_ensemble1.make_test_statistic_histogram(fit=True, model="SkewedGaussianModel")
        plt.subplot(1, 2, 2)
        hist2, fit2 = likelihood_ensemble2.make_test_statistic_histogram(fit=True, model="GaussianModel")
        plt.savefig(f"{path}/likelihood_ratios.png")
        if args.make_plots: plt.show()
        plt.close()
        
        outfile = open(f"{path}/likelihood_ratio_values.txt", "w")
        for i in range(len(hist1[0])):
            outfile.write(f"{hist1[1][i]},{hist1[0][i]},{hist2[1][i]},{hist2[0][i]}\n")
        outfile.close()
"""
