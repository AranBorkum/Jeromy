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
        
    syst = h.get_systematics(config)

    
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
            O15_s      = DataHZ.O15    .shuffle_slightly(syst["O15_syst"]    , 1, 3)
            F17_s      = DataHZ.F17    .shuffle_slightly(syst["F17_syst"]    , 1, 3)
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

    plt.subplot(1, 3, 3)
    plt.hist(test_statistic0, bins=50, label="No CNO")
    plt.hist(test_statistic1, bins=50, label="CNO")
    plt.legend(loc="best")
    ensure_directory_exists(f"{Jeromy.__jeromy_output__}/{config.systematics}_{config.metallicity}")
    outdir = f"{Jeromy.__jeromy_output__}/{config.systematics}_{config.metallicity}"
    plt.savefig(f"{outdir}exp_{config.exposure}_bkg_{config.background_reduction}_{config.metallicity}.png")
    plt.close()                 
    # plt.show()

    # plt.figure()
    # xRange = np.linspace(xValues[0], xValues[-1], 1000)
    # plt.plot(xRange, gaussian(xRange, center1, sigma1, 1))
    I_l, error = quad(gaussian, -np.inf, center2, args=(center1, sigma1, 1))
    # plt.close()
    output = f"{config.metallicity},"
    output += f"{config.background_reduction},"
    output += f"{config.exposure},"
    output += f"{2*I_l}"
    print(output)
    

