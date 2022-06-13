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

class CountingExperiment():

    def __init__(self, name, roi, config, shuffle=False, unc_table=False):
        self._name       = name
        self._low        = roi[0]
        self._high       = roi[1]
        self.config      = config
        self.shuffle     = shuffle
        self.unc_table   = unc_table
        self.syst = h.get_systematics(self.config)

        if self.config.metallicity == "high":
            if shuffle:
                self.O15      = DataHZ.O15    .shuffle_slightly(self.syst["O15_syst"]    , self._low, self._high)
                self.F17      = DataHZ.F17    .shuffle_slightly(self.syst["F17_syst"]    , self._low, self._high)
                self.B8_CC    = DataHZ.B8_CC  .shuffle_slightly(self.syst["B8_CC_syst"]  , self._low, self._high)
                self.B8_ES    = DataHZ.B8_ES  .shuffle_slightly(self.syst["B8_ES_syst"]  , self._low, self._high)
                self.HEP_CC   = DataHZ.HEP_CC .shuffle_slightly(self.syst["HEP_CC_syst"] , self._low, self._high)
                self.HEP_ES   = DataHZ.HEP_ES .shuffle_slightly(self.syst["HEP_ES_syst"] , self._low, self._high)
                self.radon    = DataHZ.Radon  .shuffle_slightly(self.syst["Radon_syst"]  , self._low, self._high)
                self.neutron  = DataHZ.Neutron.shuffle_slightly(self.syst["Neutron_syst"], self._low, self._high)
                self.argon    = DataHZ.Ar42   .shuffle_slightly(self.syst["Ar42_syst"]   , self._low, self._high)
            else:
                self.O15      = DataHZ.O15   
                self.F17      = DataHZ.F17   
                self.B8_CC    = DataHZ.B8_CC 
                self.B8_ES    = DataHZ.B8_ES 
                self.HEP_CC   = DataHZ.HEP_CC
                self.HEP_ES   = DataHZ.HEP_ES
                self.radon    = DataHZ.Radon 
                self.neutron  = DataHZ.Neutron
                self.argon    = DataHZ.Ar42  

        if self.config.metallicity == "low":
            if shuffle:
                self.O15      = DataLZ.O15    .shuffle_slightly(self.syst["O15_syst"]    , self._low, self._high)
                self.F17      = DataLZ.F17    .shuffle_slightly(self.syst["F17_syst"]    , self._low, self._high)
                self.B8_CC    = DataLZ.B8_CC  .shuffle_slightly(self.syst["B8_CC_syst"]  , self._low, self._high)
                self.B8_ES    = DataLZ.B8_ES  .shuffle_slightly(self.syst["B8_ES_syst"]  , self._low, self._high)
                self.HEP_CC   = DataLZ.HEP_CC .shuffle_slightly(self.syst["HEP_CC_syst"] , self._low, self._high)
                self.HEP_ES   = DataLZ.HEP_ES .shuffle_slightly(self.syst["HEP_ES_syst"] , self._low, self._high)
                self.radon    = DataLZ.Radon  .shuffle_slightly(self.syst["Radon_syst"]  , self._low, self._high) 
                self.neutron  = DataLZ.Neutron.shuffle_slightly(self.syst["Neutron_syst"], self._low, self._high)
                self.argon    = DataLZ.Ar42   .shuffle_slightly(self.syst["Ar42_syst"]   , self._low, self._high) 
            else:
                self.O15      = DataLZ.O15   
                self.F17      = DataLZ.F17   
                self.B8_CC    = DataLZ.B8_CC 
                self.B8_ES    = DataLZ.B8_ES 
                self.HEP_CC   = DataLZ.HEP_CC
                self.HEP_ES   = DataLZ.HEP_ES
                self.radon    = DataLZ.Radon 
                self.neutron  = DataLZ.Neutron
                self.argon    = DataLZ.Ar42  

        rate, uncertainty = self.calculate_rate()
        self.cno_rate = rate
        self.cno_uncertainty = uncertainty
            
    def _calculate_epsilon(self, spectrum, exposure, background_reduction=0):
        rate_total = 0
        rate_in_roi = 0
        epsilon = 0
        for i in spectrum.map:
            rate_total += spectrum.map[i]
            if i >= self._low and i <= self._high:
                rate_in_roi += spectrum.map[i]

        rate_total *= float(exposure); rate_in_roi *= float(exposure)
        if background_reduction:
            rate_total /= float(background_reduction)
            rate_in_roi /= float(background_reduction)
        
        epsilon = rate_in_roi / rate_total
        return rate_total, rate_in_roi, epsilon

    def calculate_rate(self, print_table=True):
        # print("Running counting experiment over an energy range "
        #       f"{self._low} - {self._high} MeV")

        O15     = self._calculate_epsilon(self.O15    , self.config.exposure)
        F17     = self._calculate_epsilon(self.F17    , self.config.exposure)
        B8_CC   = self._calculate_epsilon(self.B8_CC  , self.config.exposure)
        B8_ES   = self._calculate_epsilon(self.B8_ES  , self.config.exposure)
        HEP_CC  = self._calculate_epsilon(self.HEP_CC , self.config.exposure)
        HEP_ES  = self._calculate_epsilon(self.HEP_ES , self.config.exposure)
        radon   = self._calculate_epsilon(self.radon  , self.config.exposure, self.config.background_reduction)
        neutron = self._calculate_epsilon(self.neutron, self.config.exposure, self.config.background_reduction)
        argon   = self._calculate_epsilon(self.argon  , self.config.exposure, self.config.background_reduction)

        my_values = [O15, F17, B8_CC, B8_ES, HEP_CC, HEP_ES, radon, neutron, argon]
        total_rate_in_roi = sum(i[1] for i in my_values)
        rate = total_rate_in_roi - sum(i[0]*i[2] for i in my_values[2:])
        cno_epsilon = (O15[2] + F17[2])/2
        rate /= cno_epsilon
        uncertainty = np.sqrt(total_rate_in_roi)**2
        uncertainty += (B8_CC[0]   * B8_CC[2]   * self.syst["B8_CC_syst"]  )**2
        uncertainty += (B8_ES[0]   * B8_ES[2]   * self.syst["B8_ES_syst"]  )**2
        uncertainty += (HEP_CC[0]  * HEP_CC[2]  * self.syst["HEP_CC_syst"] )**2
        uncertainty += (HEP_ES[0]  * HEP_ES[2]  * self.syst["HEP_ES_syst"] )**2
        uncertainty += (radon[0]   * radon[2]   * self.syst["Radon_syst"]  )**2
        uncertainty += (neutron[0] * neutron[2] * self.syst["Neutron_syst"])**2
        uncertainty += (argon[0]   * argon[2]   * self.syst["Ar42_syst"]   )**2
        uncertainty = np.sqrt(uncertainty)
        uncertainty /= cno_epsilon

        if self.unc_table:
            print("B8_CC   %.2E" % (B8_CC[0]   * B8_CC[2]   * self.syst["B8_CC_syst"]  ))
            print("B8_ES   %.2E" % (B8_ES[0]   * B8_ES[2]   * self.syst["B8_ES_syst"]  ))
            print("HEP_CC  %.2E" % (HEP_CC[0]  * HEP_CC[2]  * self.syst["HEP_CC_syst"] ))
            print("HEP_ES  %.2E" % (HEP_ES[0]  * HEP_ES[2]  * self.syst["HEP_ES_syst"] ))
            print("radon   %.2E" % (radon[0]   * radon[2]   * self.syst["Radon_syst"]  ))
            print("neutron %.2E" % (neutron[0] * neutron[2] * self.syst["Neutron_syst"]))
            print("argon   %.2E" % (argon[0]   * argon[2]   * self.syst["Ar42_syst"]   ))
            


        
        
        table_B8_CC   = (B8_CC[0]   * B8_CC[2]  , np.sqrt(B8_CC[0]*B8_CC[2])    , B8_CC[0]   * B8_CC[2]   * self.syst["B8_CC_syst"]  )
        table_B8_ES   = (B8_ES[0]   * B8_ES[2]  , np.sqrt(B8_ES[0]*B8_ES[2])    , B8_ES[0]   * B8_ES[2]   * self.syst["B8_ES_syst"]  )
        table_HEP_CC  = (HEP_CC[0]  * HEP_CC[2] , np.sqrt(HEP_CC[0]*HEP_CC[2])  , HEP_CC[0]  * HEP_CC[2]  * self.syst["HEP_CC_syst"] )
        table_HEP_ES  = (HEP_ES[0]  * HEP_ES[2] , np.sqrt(HEP_ES[0]*HEP_ES[2])  , HEP_ES[0]  * HEP_ES[2]  * self.syst["HEP_ES_syst"] )
        table_radon   = (radon[0]   * radon[2]  , np.sqrt(radon[0]*radon[2])    , radon[0]   * radon[2]   * self.syst["Radon_syst"]  )
        table_neutron = (neutron[0] * neutron[2], np.sqrt(neutron[0]*neutron[2]), neutron[0] * neutron[2] * self.syst["Neutron_syst"])
        table_argon   = (argon[0]   * argon[2]  , np.sqrt(argon[0]*argon[2])    , argon[0]   * argon[2]   * self.syst["Ar42_syst"]   )


        
        if print_table:
            print("{: <15} {: <15}".format("Source", "Rate in ROI ± (stat.) ± (syst.)"))
            print("-"*45)
            # print("{: <15} {: <15}".format("CNO"     , "%.2E ± %.2E ± %.2E"% (rate, cno_stat, uncertainty)))
            print("{: <15} {: <15}".format("Boron8 cc", "%.2E ± %.2E ± %.2E"% table_B8_CC  ))
            print("{: <15} {: <15}".format("Boron8 es", "%.2E ± %.2E ± %.2E"% table_B8_ES  ))
            print("{: <15} {: <15}".format("HEP cc"   , "%.2E ± %.2E ± %.2E"% table_HEP_CC ))
            print("{: <15} {: <15}".format("HEP es"   , "%.2E ± %.2E ± %.2E"% table_HEP_ES ))
            print("{: <15} {: <15}".format("Radon"    , "%.2E ± %.2E ± %.2E"% table_radon  ))
            print("{: <15} {: <15}".format("Neutron"  , "%.2E ± %.2E ± %.2E"% table_neutron))
            print("{: <15} {: <15}".format("Argon 42" , "%.2E ± %.2E ± %.2E"% table_argon  ))
            print("-"*45)
            print("{: <15} {: <15}".format("CNO"     , "%.2E ± %.2E ± %.2E"% (rate, np.sqrt(rate), uncertainty)))
        return rate, uncertainty
        
class CountingExperimentEnsemble():

    def __init__(self, name, roi, config, ntrials=1000):
        self._name = name
        self._low  = roi[0]
        self._high = roi[1]
        self.config = config
        self.syst = h.get_systematics(self.config)
        self.ntrials = ntrials

        cno_rates = []
        cno_uncertainties = []

        print("Running Trials")
        for i in tqdm(range(self.ntrials)):
            trial = CountingExperiment("trial_%i"%i, [1, 3], self.config, shuffle=True)
            if not isnan(trial.cno_rate) and not isnan(trial.cno_uncertainty):
                cno_rates.append(trial.cno_rate)
                cno_uncertainties.append(trial.cno_uncertainty)
            
        self.rates = cno_rates
        self.uncertainties = cno_uncertainties
        
    def make_histogram(self, bins=50, fit=False):
        plt.figure(f"CNO_rate_toy_{self.config.name}")
        hist = plt.hist(self.rates, bins=bins)
        plt.xlabel(f"CNO Rate [evts / {10*float(self.config.exposure)} kt-year]")
        if fit:
            xValues = hist[1][:bins]
            yValues = hist[0][:bins]
            max_bin = np.argmax(yValues)
            model = GaussianModel()
            params = model.make_params(amplitude=yValues[max_bin], center=xValues[max_bin], sigma=1000)
            result = model.fit(yValues, params, x=xValues)
            plt.plot(xValues, result.best_fit)
        plt.show()
        plt.close()
        return hist

    def fit_histogram(self, bins=50, fit=True):
        plt.figure(f"CNO_rate_toy_{self.config.name}")
        hist = plt.hist(self.rates, bins=bins)
        plt.xlabel(f"CNO Rate [evts / {10*float(self.config.exposure)} kt-year]")
        if fit:
            xValues = hist[1][:bins]
            yValues = hist[0][:bins]
            max_bin = np.argmax(yValues)
            model = GaussianModel()
            params = model.make_params(amplitude=yValues[max_bin], center=xValues[max_bin], sigma=1000)
            result = model.fit(yValues, params, x=xValues)
            plt.plot(xValues, result.best_fit)
        plt.close()
        return result.best_values

    
    def evaluate(self):
        best_fit = self.fit_histogram()
        Spectra = true_solar.SolarTrue(self.config)
        CNO_high = Spectra.CNO_high
        CNO_low = Spectra.CNO_low

        HighZ  = f"High Z value: {CNO_high.mean} ± {CNO_high.error} "
        center = best_fit["center"]
        HighZ += f"calculated value: {center} "
        p_value_high = h.p_value_to_sigma(CNO_high.p_value(best_fit["center"]))
        HighZ += f"probability: {p_value_high}"
        print(HighZ)

        LowZ  = f"Low Z value: {CNO_low.mean} ± {CNO_low.error} "
        center = best_fit["center"]
        LowZ += f"calculated value: {center} "
        p_value_low = h.p_value_to_sigma(CNO_low.p_value(best_fit["center"]))
        LowZ += f"probability: {p_value_low}"
        print(LowZ)


class CountingExperiment_B8():

    def __init__(self, name, roi, config, shuffle=False):
        self._name       = name
        self._low        = roi[0]
        self._high       = roi[1]
        self.config      = config
        self.shuffle     = shuffle
        self.syst = h.get_systematics(self.config)

        if self.config.metallicity == "high":
            if shuffle:
                self.B8_CC    = DataHZ.B8_CC  .shuffle_slightly(self.syst["B8_CC_syst"]  , self._low, self._high)
                self.B8_ES    = DataHZ.B8_ES  .shuffle_slightly(self.syst["B8_ES_syst"]  , self._low, self._high)
                self.HEP_CC   = DataHZ.HEP_CC .shuffle_slightly(self.syst["HEP_CC_syst"] , self._low, self._high)
                self.HEP_ES   = DataHZ.HEP_ES .shuffle_slightly(self.syst["HEP_ES_syst"] , self._low, self._high)
                self.radon    = DataHZ.Radon  .shuffle_slightly(self.syst["Radon_syst"]  , self._low, self._high)
                self.neutron  = DataHZ.Neutron.shuffle_slightly(self.syst["Neutron_syst"], self._low, self._high)
                self.argon    = DataHZ.Ar42   .shuffle_slightly(self.syst["Ar42_syst"]   , self._low, self._high)
            else:
                self.B8_CC    = DataHZ.B8_CC 
                self.B8_ES    = DataHZ.B8_ES 
                self.HEP_CC   = DataHZ.HEP_CC
                self.HEP_ES   = DataHZ.HEP_ES
                self.radon    = DataHZ.Radon 
                self.neutron  = DataHZ.Neutron
                self.argon    = DataHZ.Ar42  

        if self.config.metallicity == "low":
            if shuffle:
                self.B8_CC    = DataLZ.B8_CC  .shuffle_slightly(self.syst["B8_CC_syst"]  , self._low, self._high)
                self.B8_ES    = DataLZ.B8_ES  .shuffle_slightly(self.syst["B8_ES_syst"]  , self._low, self._high)
                self.HEP_CC   = DataLZ.HEP_CC .shuffle_slightly(self.syst["HEP_CC_syst"] , self._low, self._high)
                self.HEP_ES   = DataLZ.HEP_ES .shuffle_slightly(self.syst["HEP_ES_syst"] , self._low, self._high)
                self.radon    = DataLZ.Radon  .shuffle_slightly(self.syst["Radon_syst"]  , self._low, self._high) 
                self.neutron  = DataLZ.Neutron.shuffle_slightly(self.syst["Neutron_syst"], self._low, self._high)
                self.argon    = DataLZ.Ar42   .shuffle_slightly(self.syst["Ar42_syst"]   , self._low, self._high) 
            else:
                self.B8_CC    = DataLZ.B8_CC 
                self.B8_ES    = DataLZ.B8_ES 
                self.HEP_CC   = DataLZ.HEP_CC
                self.HEP_ES   = DataLZ.HEP_ES
                self.radon    = DataLZ.Radon 
                self.neutron  = DataLZ.Neutron
                self.argon    = DataLZ.Ar42  

        rate, uncertainty = self.calculate_rate()
        self.b8_rate = rate
        self.b8_uncertainty = uncertainty
            
    def _calculate_epsilon(self, spectrum, exposure, background_reduction=0):
        rate_total = 0
        rate_in_roi = 0
        epsilon = 0
        for i in spectrum.map:
            rate_total += spectrum.map[i]
            if i >= self._low and i <= self._high:
                rate_in_roi += spectrum.map[i]

        rate_total *= float(exposure); rate_in_roi *= float(exposure)
        if background_reduction:
            rate_total /= float(background_reduction)
            rate_in_roi /= float(background_reduction)
        
        epsilon = rate_in_roi / rate_total
        return rate_total, rate_in_roi, epsilon

    def calculate_rate(self, print_table=False):
        # print("Running counting experiment over an energy range "
        #       f"{self._low} - {self._high} MeV")

        B8_CC   = self._calculate_epsilon(self.B8_CC  , self.config.exposure)
        B8_ES   = self._calculate_epsilon(self.B8_ES  , self.config.exposure)
        HEP_CC  = self._calculate_epsilon(self.HEP_CC , self.config.exposure)
        HEP_ES  = self._calculate_epsilon(self.HEP_ES , self.config.exposure)
        radon   = self._calculate_epsilon(self.radon  , self.config.exposure, self.config.background_reduction)
        neutron = self._calculate_epsilon(self.neutron, self.config.exposure, self.config.background_reduction)
        argon   = self._calculate_epsilon(self.argon  , self.config.exposure, self.config.background_reduction)

        my_values = [B8_CC, B8_ES, HEP_CC, HEP_ES, radon, neutron, argon]
        total_rate_in_roi = sum(i[1] for i in my_values)
        rate = total_rate_in_roi - sum(i[0]*i[2] for i in my_values[2:])
        b8_epsilon = (B8_CC[2] + B8_ES[2])/2
        rate /= b8_epsilon
        uncertainty = np.sqrt(total_rate_in_roi)**2
        uncertainty += (HEP_CC[0]  * HEP_CC[2]  * self.syst["HEP_CC_syst"] )**2
        uncertainty += (HEP_ES[0]  * HEP_ES[2]  * self.syst["HEP_ES_syst"] )**2
        uncertainty += (radon[0]   * radon[2]   * self.syst["Radon_syst"]  )**2
        uncertainty += (neutron[0] * neutron[2] * self.syst["Neutron_syst"])**2
        uncertainty += (argon[0]   * argon[2]   * self.syst["Ar42_syst"]   )**2
        uncertainty = np.sqrt(uncertainty)
        uncertainty /= b8_epsilon
        
        
        table_HEP_CC  = (HEP_CC[0]  * HEP_CC[2] , np.sqrt(HEP_CC[0]*HEP_CC[2])  , HEP_CC[0]  * HEP_CC[2]  * self.syst["HEP_CC_syst"] )
        table_HEP_ES  = (HEP_ES[0]  * HEP_ES[2] , np.sqrt(HEP_ES[0]*HEP_ES[2])  , HEP_ES[0]  * HEP_ES[2]  * self.syst["HEP_ES_syst"] )
        table_radon   = (radon[0]   * radon[2]  , np.sqrt(radon[0]*radon[2])    , radon[0]   * radon[2]   * self.syst["Radon_syst"]  )
        table_neutron = (neutron[0] * neutron[2], np.sqrt(neutron[0]*neutron[2]), neutron[0] * neutron[2] * self.syst["Neutron_syst"])
        table_argon   = (argon[0]   * argon[2]  , np.sqrt(argon[0]*argon[2])    , argon[0]   * argon[2]   * self.syst["Ar42_syst"]   )


        
        if print_table:
            print("{: <15} {: <15}".format("Source", "Rate in ROI ± (stat.) ± (syst.)"))
            print("-"*45)
            # print("{: <15} {: <15}".format("CNO"     , "%.2E ± %.2E ± %.2E"% (rate, cno_stat, uncertainty)))
            print("{: <15} {: <15}".format("HEP cc"   , "%.2E ± %.2E ± %.2E"% table_HEP_CC ))
            print("{: <15} {: <15}".format("HEP es"   , "%.2E ± %.2E ± %.2E"% table_HEP_ES ))
            print("{: <15} {: <15}".format("Radon"    , "%.2E ± %.2E ± %.2E"% table_radon  ))
            print("{: <15} {: <15}".format("Neutron"  , "%.2E ± %.2E ± %.2E"% table_neutron))
            print("{: <15} {: <15}".format("Argon 42" , "%.2E ± %.2E ± %.2E"% table_argon  ))
            print("-"*45)
            print("{: <15} {: <15}".format("B8"     , "%.2E ± %.2E ± %.2E"% (rate, np.sqrt(rate), uncertainty)))
        return rate, uncertainty
