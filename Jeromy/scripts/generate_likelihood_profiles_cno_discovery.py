#!/usr/bin/env python3

import os
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
import Jeromy.analyzers.cno_discovery        as cno_discovery

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
   

if  __name__ == "__main__":

    os.system("clear")
    args = create_parser()
    if args.configuration:
        config = configuration_loader.ConfigurationLoader(args.configuration)
    else:
        config = configuration_loader.ConfigurationLoader(systematics=args.systematics,
                                                          metallicity=args.metallicity,
                                                          exposure=args.exposure,
                                                          background_reduction=args.background_reduction)
        
    syst = h.get_systematics(config)
    CNO_Discovery = cno_discovery.DiscoveryCNO(config, [1, 3], 1000)
    
    bottom_value = min(min(CNO_Discovery.cno_injected_likelihood[0]),
                       min(CNO_Discovery.background_only_likelihood[0]))
    top_value    = max(max(CNO_Discovery.cno_injected_likelihood[0]),
                       max(CNO_Discovery.background_only_likelihood[0]))

    q_value_range = np.linspace(bottom_value, top_value, 10000)
    
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.hist(CNO_Discovery.background.xValues, bins=range(9),
             weights=CNO_Discovery.background.yValues, histtype="step",
             label="Backgrounds only")
    plt.hist(CNO_Discovery.signal.xValues    , bins=range(9),
             weights=CNO_Discovery.signal.yValues, histtype="step",
             label="CNO + Backgrounds")
    plt.yscale("log")
    plt.ylim([10, 10**9])
    plt.xlabel("Reconstructed neutrino energy [MeV]", fontsize=16)
    plt.ylabel("Rate [evts / 10 kt-year / MeV]", fontsize=16)
    plt.legend(loc="best")
    
    plt.subplot(1, 2, 2)
    # plt.hist(CNO_Discovery.cno_injected_likelihood[0], bins=CNO_Discovery.cno_injected_likelihood[0],
    #          weights=CNO_Discovery.cno_injected_likelihood[1],
    #          label="Signal + Background", histtype="step")
    plt.hist(CNO_Discovery.background_only_likelihood[0], bins=CNO_Discovery.background_only_likelihood[0],
             weights=CNO_Discovery.background_only_likelihood[1],
             label="Background", histtype="step")
    plt.plot(CNO_Discovery.background_only_likelihood[0], CNO_Discovery.background_fit.best_fit)
    # plt.plot(CNO_Discovery.cno_injected_likelihood[0], CNO_Discovery.signal_fit.best_fit)

    # plt.legend(loc="best")
    plt.xlabel("Test Statistic", fontsize=16)
    title =  f"exposure: {float(config.exposure) * 10}, bkg reduction: {config.background_reduction}x "
    title += f"Z: {config.metallicity}"
    plt.title(title)


    # Find the point where the background likelihood reaches 0.003
    background_fit_center = CNO_Discovery.background_fit.best_values["center"]
    background_fit_sigma  = CNO_Discovery.background_fit.best_values["sigma"]
    
    q_value_required_1sigma = background_fit_center + (1 * background_fit_sigma)
    q_value_required_2sigma = background_fit_center + (2 * background_fit_sigma)
    q_value_required_3sigma = background_fit_center + (3 * background_fit_sigma)
    q_value_required_4sigma = background_fit_center + (4 * background_fit_sigma)
    q_value_required_5sigma = background_fit_center + (5 * background_fit_sigma)
    q_value_required_6sigma = background_fit_center + (6 * background_fit_sigma)


    

    total_backgrounds = 0
    total_signal      = 0
    for i in CNO_Discovery.background.map:
        if i >= CNO_Discovery.roi[0] and i <= CNO_Discovery.roi[1]:
            total_backgrounds += CNO_Discovery.background.map[i]
            total_signal      += CNO_Discovery.signal.map[i]

    value = 1.40E+05
    unc = 0.52
    q_estimate1 = h.poisson_likelihood(total_backgrounds, total_backgrounds+(value*(1-unc)))
    q_estimate2 = h.poisson_likelihood(total_backgrounds, total_backgrounds+(value*(1    )))
    q_estimate3 = h.poisson_likelihood(total_backgrounds, total_backgrounds+(value*(1+unc)))
    plt.plot([q_estimate1, q_estimate1], [0, 0.06], "k--")
    plt.plot([q_estimate2, q_estimate2], [0, 0.06], "k--")
    plt.plot([q_estimate3, q_estimate3], [0, 0.06], "k--")
    I_1 = quad(h.gaussian, q_estimate1, np.inf, args=(background_fit_center, background_fit_sigma, 1))
    I_2 = quad(h.gaussian, q_estimate2, np.inf, args=(background_fit_center, background_fit_sigma, 1))
    I_3 = quad(h.gaussian, q_estimate3, np.inf, args=(background_fit_center, background_fit_sigma, 1))

    print(erfinv(1-I_1[0]) * np.sqrt(2))
    print(erfinv(1-I_2[0]) * np.sqrt(2))
    print(erfinv(1-I_3[0]) * np.sqrt(2))
    
            
    yValues = np.linspace(total_backgrounds, total_signal*2, 1000)
    xValues = [h.poisson_likelihood(total_backgrounds, i) for i in yValues]
    yValues = [i-total_backgrounds for i in yValues]

    n_events_required_1sigma = 0
    n_events_required_2sigma = 0
    n_events_required_3sigma = 0
    n_events_required_4sigma = 0
    n_events_required_5sigma = 0
    n_events_required_6sigma = 0
            
    for i, x in enumerate(xValues):
        if x > q_value_required_1sigma and not n_events_required_1sigma:
            n_events_required_1sigma = yValues[i]
        if x > q_value_required_2sigma and not n_events_required_2sigma:
            n_events_required_2sigma = yValues[i]
        if x > q_value_required_3sigma and not n_events_required_3sigma:
            n_events_required_3sigma = yValues[i]
        if x > q_value_required_4sigma and not n_events_required_4sigma:
            n_events_required_4sigma = yValues[i]
        if x > q_value_required_5sigma and not n_events_required_5sigma:
            n_events_required_5sigma = yValues[i]
        if x > q_value_required_6sigma and not n_events_required_6sigma:
            n_events_required_6sigma = yValues[i]

        if (n_events_required_2sigma and n_events_required_3sigma and n_events_required_4sigma and
            n_events_required_5sigma and n_events_required_6sigma):
            break

    
    print(n_events_required_1sigma,
          n_events_required_2sigma,
          n_events_required_3sigma,
          n_events_required_4sigma)
    # plt.subplot(1, 3, 3)    
    # plt.plot(xValues, yValues)

    # plt.fill_between([0, q_value_required_6sigma],
    #                  [n_events_required_6sigma, n_events_required_6sigma],
    #                  [0, 0], alpha=0.2)

    # plt.fill_between([0, q_value_required_5sigma],
    #                  [n_events_required_5sigma, n_events_required_5sigma],
    #                  [0, 0], alpha=0.2)

    # plt.fill_between([0, q_value_required_4sigma],
    #                  [n_events_required_4sigma, n_events_required_4sigma],
    #                  [0, 0], alpha=0.2)


    # plt.fill_between([0, q_value_required_3sigma],
    #                  [n_events_required_3sigma, n_events_required_3sigma],
    #                  [0, 0], alpha=0.2)

    # plt.fill_between([0, q_value_required_2sigma],
    #                  [n_events_required_2sigma, n_events_required_2sigma],
    #                  [0, 0], alpha=0.2)

    # plt.fill_between([0, q_value_required_1sigma],
    #                  [n_events_required_1sigma, n_events_required_1sigma],
    #                  [0, 0], alpha=0.0)

    
    # plt.plot([0, q_value_required_1sigma], [n_events_required_1sigma, n_events_required_1sigma], "k--")
    # plt.plot([q_value_required_1sigma, q_value_required_1sigma], [0, n_events_required_1sigma], "k--")
    
    # plt.plot([0, q_value_required_2sigma], [n_events_required_2sigma, n_events_required_2sigma], "k--")
    # plt.plot([q_value_required_2sigma, q_value_required_2sigma], [0, n_events_required_2sigma], "k--")

    # plt.plot([0, q_value_required_3sigma], [n_events_required_3sigma, n_events_required_3sigma], "k--")
    # plt.plot([q_value_required_3sigma, q_value_required_3sigma], [0, n_events_required_3sigma], "k--")

    # plt.plot([0, q_value_required_4sigma], [n_events_required_4sigma, n_events_required_4sigma], "k--")
    # plt.plot([q_value_required_4sigma, q_value_required_4sigma], [0, n_events_required_4sigma], "k--")

    # plt.plot([0, q_value_required_5sigma], [n_events_required_5sigma, n_events_required_5sigma], "k--")
    # plt.plot([q_value_required_5sigma, q_value_required_5sigma], [0, n_events_required_5sigma], "k--")

    # plt.plot([0, q_value_required_6sigma], [n_events_required_6sigma, n_events_required_6sigma], "k--")
    # plt.plot([q_value_required_6sigma, q_value_required_6sigma], [0, n_events_required_6sigma], "k--")

    # plt.xlabel("test statistic", fontsize=16)
    # plt.ylabel("number of events required", fontsize=16)
    # plt.title("")
    # plt.ylim([0, n_events_required_5sigma*2])
    # plt.xlim([0, q_value_required_5sigma*2])
    # plt.xticks(rotation=30)
    # # plt.close()

    # plt.errorbar(5000, 2.5E+07, yerr=2.5E+07*0.62, fmt="x")
    plt.show()



    
    """
    if config.metallicity == "high":
        O15      = DataHZ.O15
        F17      = DataHZ.F17
        B8_CC    = DataHZ.B8_CC 
        B8_ES    = DataHZ.B8_ES 
        HEP_CC   = DataHZ.HEP_CC
        HEP_ES   = DataHZ.HEP_ES
        radon    = DataHZ.Radon 
        neutron  = DataHZ.Neutron
        argon    = DataHZ.Ar42  
        
    if config.metallicity == "low":
        O15      = DataLZ.O15
        F17      = DataLZ.F17
        B8_CC    = DataLZ.B8_CC 
        B8_ES    = DataLZ.B8_ES 
        HEP_CC   = DataLZ.HEP_CC
        HEP_ES   = DataLZ.HEP_ES
        radon    = DataLZ.Radon 
        neutron  = DataLZ.Neutron
        argon    = DataLZ.Ar42  

    scale_spectrum(O15    , config.exposure, 1)
    scale_spectrum(F17    , config.exposure, 1)    
    scale_spectrum(B8_CC  , config.exposure, 1)
    scale_spectrum(B8_ES  , config.exposure, 1)    
    scale_spectrum(HEP_CC , config.exposure, 1)
    scale_spectrum(HEP_ES , config.exposure, 1)
    scale_spectrum(radon  , config.exposure, config.background_reduction)
    scale_spectrum(neutron, config.exposure, config.background_reduction)
    scale_spectrum(argon  , config.exposure, config.background_reduction)
        
    total = FI.ImportFile("TotalBackground")
    total.combine_inputs([radon, neutron, argon, B8_CC, B8_ES, HEP_CC, HEP_ES])
    baseline = counting_experiment.CountingExperiment("baseline", [1, 3], config)

    test_statistic0 = []
    test_statistic1 = []
    for i in tqdm(range(1000)):
        if config.metallicity == "high":
            O15_s      = DataHZ.O15    .shuffle_slightly(syst["O15_syst"]    , 1, 3)
            F17_s      = DataHZ.F17    .shuffle_slightly(syst["F17_syst"]    , 1, 3)
            B8_CC_s    = DataHZ.B8_CC  .shuffle_slightly(syst["B8_CC_syst"]  , 1, 3)
            B8_ES_s    = DataHZ.B8_ES  .shuffle_slightly(syst["B8_ES_syst"]  , 1, 3)
            HEP_CC_s   = DataHZ.HEP_CC .shuffle_slightly(syst["HEP_CC_syst"] , 1, 3)
            HEP_ES_s   = DataHZ.HEP_ES .shuffle_slightly(syst["HEP_ES_syst"] , 1, 3)
            radon_s    = DataHZ.Radon  .shuffle_slightly(syst["Radon_syst"]  , 1, 3)
            neutron_s  = DataHZ.Neutron.shuffle_slightly(syst["Neutron_syst"], 1, 3)
            argon_s    = DataHZ.Ar42   .shuffle_slightly(syst["Ar42_syst"]   , 1, 3)
        if config.metallicity == "low":
            O15_s      = DataLZ.O15    .shuffle_slightly(syst["O15_syst"]    , 1, 3)
            F17_s      = DataLZ.F17    .shuffle_slightly(syst["F17_syst"]    , 1, 3)
            B8_CC_s    = DataLZ.B8_CC  .shuffle_slightly(syst["B8_CC_syst"]  , 1, 3)
            B8_ES_s    = DataLZ.B8_ES  .shuffle_slightly(syst["B8_ES_syst"]  , 1, 3)
            HEP_CC_s   = DataLZ.HEP_CC .shuffle_slightly(syst["HEP_CC_syst"] , 1, 3)
            HEP_ES_s   = DataLZ.HEP_ES .shuffle_slightly(syst["HEP_ES_syst"] , 1, 3)
            radon_s    = DataLZ.Radon  .shuffle_slightly(syst["Radon_syst"]  , 1, 3)
            neutron_s  = DataLZ.Neutron.shuffle_slightly(syst["Neutron_syst"], 1, 3)
            argon_s    = DataLZ.Ar42   .shuffle_slightly(syst["Ar42_syst"]   , 1, 3)

        scale_spectrum(O15_s    , config.exposure, 1)
        scale_spectrum(F17_s    , config.exposure, 1)    
        scale_spectrum(B8_CC_s  , config.exposure, 1)
        scale_spectrum(B8_ES_s  , config.exposure, 1)    
        scale_spectrum(HEP_CC_s , config.exposure, 1)
        scale_spectrum(HEP_ES_s , config.exposure, 1)
        scale_spectrum(radon_s  , config.exposure, config.background_reduction)
        scale_spectrum(neutron_s, config.exposure, config.background_reduction)
        scale_spectrum(argon_s  , config.exposure, config.background_reduction)
        total_s0 = FI.ImportFile("TotalBackground_s")
        total_s1 = FI.ImportFile("TotalBackground_s")
        total_s0.combine_inputs([radon_s, neutron_s, argon_s, B8_CC_s, B8_ES_s, HEP_CC_s, HEP_ES_s])
        total_s1.combine_inputs([radon_s, neutron_s, argon_s, B8_CC_s,
                                 B8_ES_s, HEP_CC_s, HEP_ES_s, O15_s, F17_s])

        total_scaler = sum(total.yValues)


        q = 0
        for i in total.map:
            if i in total_s0.map:
                y = total.map[i] / total_scaler
                n = total_s0.map[i] / total_scaler
                q += y-n+(n*np.log(n/y))
        test_statistic0.append(q*2)

        q = 0
        for i in total.map:
            if i in total_s1.map:
                y = total.map[i] / total_scaler
                n = total_s1.map[i] / total_scaler
                q += y-n+(n*np.log(n/y))
        test_statistic1.append(q*2)

    plt.figure()
    plt.subplot(1, 3, 1)
    hist = plt.hist(test_statistic0, bins=50)
    xValues = hist[1][:-1]
    yValues = hist[0]
    model = GaussianModel()
    result = model.fit(yValues, x=xValues)
    center1 = result.best_values["center"]
    sigma1  = result.best_values["sigma" ]
    print(center1, sigma1)
    plt.plot(xValues, result.best_fit)
    
    plt.title(r"$\log(\tilde{H}_{0}/H_{0})$")
    plt.subplot(1, 3, 2)
    plt.title(r"$\log(H_{1}/H_{0})$")
    hist = plt.hist(test_statistic1, bins=50)
    xValues1 = hist[1][:-1]
    yValues1 = hist[0]
    model = GaussianModel()
    result = model.fit(yValues1, x=xValues1)
    center2 = result.best_values["center"]
    plt.plot(xValues1, result.best_fit)
    print(center2)

    plt.subplot(1, 3, 3)
    plt.hist(test_statistic0, bins=50, label="No CNO")
    plt.hist(test_statistic1, bins=50, label="CNO")
    plt.legend(loc="best")
    h.ensure_directory_exists(f"{Jeromy.__jeromy_output__}/{args.systematics}_{config.metallicity}")
    outdir = f"{Jeromy.__jeromy_output__}/{args.systematics}_{config.metallicity}"
    print(f"{outdir}/exp_{config.exposure}_bkg_{config.background_reduction}_{config.metallicity}.png")
    plt.savefig(f"{outdir}/exp_{config.exposure}_bkg_{config.background_reduction}_{config.metallicity}.png")
    plt.close()                 
    # plt.show()

    # plt.figure()
    # xRange = np.linspace(xValues[0], xValues[-1], 1000)
    # plt.plot(xRange, gaussian(xRange, center1, sigma1, 1))
    I_l, error = quad(gaussian, -np.inf, center2, args=(center1, sigma1, 1))



    
    # plt.close()
    
    # output =  f"{config.metallicity},"
    # output += f"{config.background_reduction},"
    # output += f"{config.exposure},"
    # output += f"{2*I_l},"
    # output += f"{baseline.cno_rate},"
    # output += f"{baseline.cno_uncertainty},"
    # output += f"{int(100*baseline.cno_uncertainty/baseline.cno_rate)}"
    # print(output)
    
    """
