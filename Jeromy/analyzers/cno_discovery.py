import numpy as np
import matplotlib.pyplot as plt

import Jeromy.core.true_solar as true_solar
import Jeromy.IO.FileImporter as FI
import Jeromy.IO.ImportDataHighZ as DataHZ
import Jeromy.IO.ImportDataLowZ as DataLZ
import Jeromy.utilities.helpers as h

from tqdm import tqdm
from math import isnan
from lmfit.models import GaussianModel

class DiscoveryCNO():

    def __init__(self, config, roi, nToys=10000):
        self.config = config
        self.roi    = roi
        self.nToys  = nToys
        self.syst   = h.get_systematics(config)

        if self.config.metallicity == "high":
            self.O15      = DataHZ.O15     .scale(self.config.exposure, 1)                               
            self.F17      = DataHZ.F17     .scale(self.config.exposure, 1)                               
            self.B8_CC    = DataHZ.B8_CC   .scale(self.config.exposure, 1)                               
            self.B8_ES    = DataHZ.B8_ES   .scale(self.config.exposure, 1)                               
            self.HEP_CC   = DataHZ.HEP_CC  .scale(self.config.exposure, 1)                               
            self.HEP_ES   = DataHZ.HEP_ES  .scale(self.config.exposure, 1)                               
            self.radon    = DataHZ.Radon   .scale(self.config.exposure, self.config.background_reduction)
            self.neutron  = DataHZ.Neutron .scale(self.config.exposure, self.config.background_reduction)
            self.argon    = DataHZ.Ar42    .scale(self.config.exposure, self.config.background_reduction)
        if self.config.metallicity == "low":
            self.O15      = DataLZ.O15     .scale(self.config.exposure, 1)                               
            self.F17      = DataLZ.F17     .scale(self.config.exposure, 1)                               
            self.B8_CC    = DataLZ.B8_CC   .scale(self.config.exposure, 1)                               
            self.B8_ES    = DataLZ.B8_ES   .scale(self.config.exposure, 1)                               
            self.HEP_CC   = DataLZ.HEP_CC  .scale(self.config.exposure, 1)                               
            self.HEP_ES   = DataLZ.HEP_ES  .scale(self.config.exposure, 1)                               
            self.radon    = DataLZ.Radon   .scale(self.config.exposure, self.config.background_reduction)
            self.neutron  = DataLZ.Neutron .scale(self.config.exposure, self.config.background_reduction)
            self.argon    = DataLZ.Ar42    .scale(self.config.exposure, self.config.background_reduction)

        self.background = FI.ImportFile(name="Background")
        self.signal     = FI.ImportFile(name="Signal"    )
            
        self.background.combine_inputs([self.B8_CC,
                                        self.B8_ES,
                                        self.HEP_CC,
                                        self.HEP_ES,
                                        self.radon,
                                        self.neutron,
                                        self.argon])
        self.signal.combine_inputs([self.B8_CC,
                                    self.B8_ES,
                                    self.HEP_CC,
                                    self.HEP_ES,
                                    self.radon,
                                    self.neutron,
                                    self.argon,
                                    self.O15,
                                    self.F17])


        cno_injected_likelihood    = self.run_toys()
        background_only_likelihood = self.run_toys(False)
        cno_hist = plt.hist(cno_injected_likelihood, bins=50)
        bkg_hist = plt.hist(background_only_likelihood, bins=50)
        plt.close()

        cno_normalised_likelihoods = [i/sum(cno_hist[0]) for i in cno_hist[0]]
        cno_bin_values             = cno_hist[1][:-1]
        bkg_normalised_likelihoods = [i/sum(bkg_hist[0]) for i in bkg_hist[0]]
        bkg_bin_values             = bkg_hist[1][:-1]

        self.cno_injected_likelihood = [cno_bin_values, cno_normalised_likelihoods]
        self.background_only_likelihood = [bkg_bin_values, bkg_normalised_likelihoods]
        self.background_fit = self.fit_model()[0]
        self.signal_fit     = self.fit_model()[1]


    def mean_and_sigma(self, x_values, y_values):
        max_value = 0
        mean      = 0
        sigma     = 0

        for i, x in enumerate(y_values):
            if x > max_value:
                max_value = x
                mean      = x_values[i]

        for i in y_values: sigma += (i - mean)**2
        sigma -= len(y_values)
        sigma = np.sqrt(sigma)

        return mean, sigma


                

    def fit_model(self):

        model_b = GaussianModel()
        mean_b, sigma_b = self.mean_and_sigma(*self.background_only_likelihood)
        param_b = model_b.make_params(center=mean_b, sigma=sigma_b, amplitude=1)
        result_b = model_b.fit(self.background_only_likelihood[1], param_b, x=self.background_only_likelihood[0])

        model_s = GaussianModel()
        mean_s, sigma_s = self.mean_and_sigma(*self.cno_injected_likelihood)
        param_s = model_s.make_params(center=mean_s, sigma=sigma_s, amplitude=1)
        result_s = model_s.fit(self.cno_injected_likelihood[1], param_s, x=self.cno_injected_likelihood[0])

        return result_b, result_s

        

    def poisson_likelihood(self, signal, background):

        tot_signal     = 0
        tot_background = 0
        for i in signal.map:
            if i >= self.roi[0] and i <= self.roi[1]:
                tot_signal += signal.map[i]

        for i in background.map:
            if i >= self.roi[0] and i <= self.roi[1]:
                tot_background += background.map[i]


        return 2 * (tot_background - tot_signal + tot_signal*np.log(tot_signal/tot_background))
                
        
        # test_stat = 0
        # for i in signal.map:
        #     if i >= self.roi[0] and i <= self.roi[1] and i in background.map and signal.map[i]:
        #         expected   = signal.map[i]
        #         observed   = background.map[i]
        #         _test_stat = expected - observed + (observed*np.log(observed/expected))
        #         test_stat += _test_stat
        # return test_stat * 2

                
    def run_toys(self, SignalToys=True):

        test_statistic_values = []
        if SignalToys:
            print("Generating likelihood profile for CNO injected spectrum")
            for i in range(self.nToys):
                temp_O15     = self.O15    .shuffle_slightly(self.syst["O15_syst"    ], self.roi[0], self.roi[1])
                temp_F17     = self.F17    .shuffle_slightly(self.syst["F17_syst"    ], self.roi[0], self.roi[1])
                temp_B8_CC   = self.B8_CC  .shuffle_slightly(self.syst["B8_CC_syst"  ], self.roi[0], self.roi[1])
                temp_B8_ES   = self.B8_ES  .shuffle_slightly(self.syst["B8_ES_syst"  ], self.roi[0], self.roi[1])
                temp_HEP_CC  = self.HEP_CC .shuffle_slightly(self.syst["HEP_CC_syst" ], self.roi[0], self.roi[1])
                temp_HEP_ES  = self.HEP_ES .shuffle_slightly(self.syst["HEP_ES_syst" ], self.roi[0], self.roi[1])
                temp_radon   = self.radon  .shuffle_slightly(self.syst["Radon_syst"  ], self.roi[0], self.roi[1])
                temp_neutron = self.neutron.shuffle_slightly(self.syst["Neutron_syst"], self.roi[0], self.roi[1])
                temp_argon   = self.argon  .shuffle_slightly(self.syst["Ar42_syst"   ], self.roi[0], self.roi[1])

                temp_signal = FI.ImportFile(name="temp_signal")
                temp_signal.combine_inputs([temp_O15    ,
                                            temp_F17    ,
                                            temp_B8_CC  ,
                                            temp_B8_ES  ,
                                            temp_HEP_CC ,
                                            temp_HEP_ES ,
                                            temp_radon  ,
                                            temp_neutron,
                                            temp_argon  ])
                test_statistic_values.append(self.poisson_likelihood(temp_signal, self.background))


        else:
            print("Generating likelihood profile for background only spectrum")            
            for i in range(self.nToys):
                temp_B8_CC   = self.B8_CC  .shuffle_slightly(self.syst["B8_CC_syst"  ], self.roi[0], self.roi[1])
                temp_B8_ES   = self.B8_ES  .shuffle_slightly(self.syst["B8_ES_syst"  ], self.roi[0], self.roi[1])
                temp_HEP_CC  = self.HEP_CC .shuffle_slightly(self.syst["HEP_CC_syst" ], self.roi[0], self.roi[1])
                temp_HEP_ES  = self.HEP_ES .shuffle_slightly(self.syst["HEP_ES_syst" ], self.roi[0], self.roi[1])
                temp_radon   = self.radon  .shuffle_slightly(self.syst["Radon_syst"  ], self.roi[0], self.roi[1])
                temp_neutron = self.neutron.shuffle_slightly(self.syst["Neutron_syst"], self.roi[0], self.roi[1])
                temp_argon   = self.argon  .shuffle_slightly(self.syst["Ar42_syst"   ], self.roi[0], self.roi[1])

                temp_signal = FI.ImportFile(name="temp_signal")
                temp_signal.combine_inputs([temp_B8_CC  ,
                                            temp_B8_ES  ,
                                            temp_HEP_CC ,
                                            temp_HEP_ES ,
                                            temp_radon  ,
                                            temp_neutron,
                                            temp_argon  ])
                test_statistic_values.append(self.poisson_likelihood(temp_signal, self.background))
        return test_statistic_values
