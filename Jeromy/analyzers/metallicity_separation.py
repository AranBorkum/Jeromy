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

class MetallicitySplit():

    def __init__(self, config, roi, nToys=1000):
        self.config       = config
        self.metallicity  = self.config.metallicity
        self.nToys        = nToys
        self.roi          = roi
        self.syst         = h.get_systematics(config)
        
        if self.metallicity == "high":
            self.metallicity2 = "low"
            self.O15      = DataHZ.O15     .scale(self.config.exposure, 1)                               
            self.F17      = DataHZ.F17     .scale(self.config.exposure, 1)                               
            self.B8_CC    = DataHZ.B8_CC   .scale(self.config.exposure, 1)                               
            self.B8_ES    = DataHZ.B8_ES   .scale(self.config.exposure, 1)                               
            self.HEP_CC   = DataHZ.HEP_CC  .scale(self.config.exposure, 1)                               
            self.HEP_ES   = DataHZ.HEP_ES  .scale(self.config.exposure, 1)                               
            self.radon    = DataHZ.Radon   .scale(self.config.exposure, self.config.background_reduction)
            self.neutron  = DataHZ.Neutron .scale(self.config.exposure, self.config.background_reduction)
            self.argon    = DataHZ.Ar42    .scale(self.config.exposure, self.config.background_reduction)

            self.h1_O15      = DataLZ.O15     .scale(self.config.exposure, 1)                               
            self.h1_F17      = DataLZ.F17     .scale(self.config.exposure, 1)                               
            self.h1_B8_CC    = DataLZ.B8_CC   .scale(self.config.exposure, 1)                               
            self.h1_B8_ES    = DataLZ.B8_ES   .scale(self.config.exposure, 1)                               
            self.h1_HEP_CC   = DataLZ.HEP_CC  .scale(self.config.exposure, 1)                               
            self.h1_HEP_ES   = DataLZ.HEP_ES  .scale(self.config.exposure, 1)                               
            self.h1_radon    = DataLZ.Radon   .scale(self.config.exposure, self.config.background_reduction)
            self.h1_neutron  = DataLZ.Neutron .scale(self.config.exposure, self.config.background_reduction)
            self.h1_argon    = DataLZ.Ar42    .scale(self.config.exposure, self.config.background_reduction)

            
        else:
            self.metallicity2 = "high"
            self.O15      = DataLZ.O15     .scale(self.config.exposure, 1)                               
            self.F17      = DataLZ.F17     .scale(self.config.exposure, 1)                               
            self.B8_CC    = DataLZ.B8_CC   .scale(self.config.exposure, 1)                               
            self.B8_ES    = DataLZ.B8_ES   .scale(self.config.exposure, 1)                               
            self.HEP_CC   = DataLZ.HEP_CC  .scale(self.config.exposure, 1)                               
            self.HEP_ES   = DataLZ.HEP_ES  .scale(self.config.exposure, 1)                               
            self.radon    = DataLZ.Radon   .scale(self.config.exposure, self.config.background_reduction)
            self.neutron  = DataLZ.Neutron .scale(self.config.exposure, self.config.background_reduction)
            self.argon    = DataLZ.Ar42    .scale(self.config.exposure, self.config.background_reduction)

            self.h1_O15      = DataHZ.O15     .scale(self.config.exposure, 1)                               
            self.h1_F17      = DataHZ.F17     .scale(self.config.exposure, 1)                               
            self.h1_B8_CC    = DataHZ.B8_CC   .scale(self.config.exposure, 1)                               
            self.h1_B8_ES    = DataHZ.B8_ES   .scale(self.config.exposure, 1)                               
            self.h1_HEP_CC   = DataHZ.HEP_CC  .scale(self.config.exposure, 1)                               
            self.h1_HEP_ES   = DataHZ.HEP_ES  .scale(self.config.exposure, 1)                               
            self.h1_radon    = DataHZ.Radon   .scale(self.config.exposure, self.config.background_reduction)
            self.h1_neutron  = DataHZ.Neutron .scale(self.config.exposure, self.config.background_reduction)
            self.h1_argon    = DataHZ.Ar42    .scale(self.config.exposure, self.config.background_reduction)
            

        self.signal      = FI.ImportFile(name="Signal"     )
        self.signal.combine_inputs([self.O15    ,
                                    self.F17    ,
                                    self.B8_CC  ,
                                    self.B8_ES  ,
                                    self.HEP_CC ,
                                    self.HEP_ES ,
                                    self.radon  ,
                                    self.neutron,
                                    self.argon  ])

        total_signal = 0
        for i in self.signal.map:
            if i >= self.roi[0] and i <= self.roi[1]:
                total_signal += self.signal.map[i]

        self.total_signal = total_signal
        null_hypothesis = self.run_toys()
        test_hypothesis = self.run_toys(False)
        null_hist = plt.hist(null_hypothesis, bins=50)
        test_hist = plt.hist(test_hypothesis, bins=50)
        self.null_hist = null_hist
        self.test_hist = test_hist
        plt.close()

        null_normalised_likelihoods = [i/sum(null_hist[0]) for i in null_hist[0]]
        null_bin_values             = null_hist[1][:-1]
        test_normalised_likelihoods = [i/sum(test_hist[0]) for i in test_hist[0]]
        test_bin_values             = test_hist[1][:-1]




        
        self.null_likelihoods = [null_bin_values, null_normalised_likelihoods]
        self.test_likelihoods = [test_bin_values, test_normalised_likelihoods]
        self.null_fit = self.fit_model()[0]
        self.test_fit = self.fit_model()[1]
        


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
        mean_b, sigma_b = self.mean_and_sigma(*self.null_likelihoods)
        param_b = model_b.make_params(center=mean_b, sigma=sigma_b, amplitude=1)
        result_b = model_b.fit(self.null_likelihoods[1], param_b, x=self.null_likelihoods[0])
        # result_b = model_b.fit(self.null_likelihoods[1], x=self.null_likelihoods[0])

        model_s = GaussianModel()
        mean_s, sigma_s = self.mean_and_sigma(*self.test_likelihoods)
        param_s = model_s.make_params(center=mean_s, sigma=sigma_s, amplitude=1)
        result_s = model_s.fit(self.test_likelihoods[1], param_s, x=self.test_likelihoods[0])
        # result_s = model_s.fit(self.test_likelihoods[1], x=self.test_likelihoods[0])
        
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

        
    def run_toys(self, H0=True):

        test_statistic_values = []
        if H0:
            print("Generating likelihood profile for null hypothesis")
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

                test_statistic_values.append(self.poisson_likelihood(temp_signal, self.signal))

        if not H0:
            print("Generating likelihood profile for test hypothesis")
            for i in range(self.nToys):
                temp_O15     = self.h1_O15    .shuffle_slightly(self.syst["O15_syst"    ], self.roi[0], self.roi[1])
                temp_F17     = self.h1_F17    .shuffle_slightly(self.syst["F17_syst"    ], self.roi[0], self.roi[1])
                temp_B8_CC   = self.h1_B8_CC  .shuffle_slightly(self.syst["B8_CC_syst"  ], self.roi[0], self.roi[1])
                temp_B8_ES   = self.h1_B8_ES  .shuffle_slightly(self.syst["B8_ES_syst"  ], self.roi[0], self.roi[1])
                temp_HEP_CC  = self.h1_HEP_CC .shuffle_slightly(self.syst["HEP_CC_syst" ], self.roi[0], self.roi[1])
                temp_HEP_ES  = self.h1_HEP_ES .shuffle_slightly(self.syst["HEP_ES_syst" ], self.roi[0], self.roi[1])
                temp_radon   = self.h1_radon  .shuffle_slightly(self.syst["Radon_syst"  ], self.roi[0], self.roi[1])
                temp_neutron = self.h1_neutron.shuffle_slightly(self.syst["Neutron_syst"], self.roi[0], self.roi[1])
                temp_argon   = self.h1_argon  .shuffle_slightly(self.syst["Ar42_syst"   ], self.roi[0], self.roi[1])

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
                test_statistic_values.append(self.poisson_likelihood(temp_signal, self.signal))
        return test_statistic_values
