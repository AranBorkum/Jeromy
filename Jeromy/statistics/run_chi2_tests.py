import random
import numpy as np

import Jeromy.statistics.test_statistic as ts
import Jeromy.IO.ImportDataHighZ as IDh
import Jeromy.IO.ImportDataLowZ as IDl
import Jeromy.utilities.helpers as h

from tqdm import tqdm

def shuffle_slightly(map, error, low, high):
    output = {}
    for i in map:
        if i >= low and i <= high:
            stat_error = np.sqrt(map[i])
            syst_error = map[i] * error
            total_error = np.sqrt(stat_error*stat_error + syst_error*syst_error)
            output[i] = map[i] + ((stat_error + syst_error) * random.randint(-1, 2))
    return output

class Chi2Test():

    def __init__(self, configuration, low, high, name, test=ts.BakerCousinsChi(per_bin=True), ntests=10000):
        self.test      = test
        self.test_name = name
        self.ntests    = ntests
        self.config    = configuration
        self.syst_unc  = configuration.systematic_uncertainties
        self.low       = low
        self.high      = high
        
    def run_tests(self):
        chi2_values = []
        O15_tot_syst_unc     = [self.syst_unc["O15_flux_e"], self.syst_unc["Sur_prob_e"],
                                self.syst_unc["trigger_e"], self.syst_unc["reconstruction"]]
        F17_tot_syst_unc     = [self.syst_unc["F17_flux_e"], self.syst_unc["Sur_prob_e"],
                                self.syst_unc["trigger_e"], self.syst_unc["reconstruction"]]
        B8_CC_tot_syst_unc   = [self.syst_unc["B8_flux_e"], self.syst_unc["Sur_prob_e"],
                                self.syst_unc["trigger_e"], self.syst_unc["reconstruction"]]
        B8_ES_tot_syst_unc   = [self.syst_unc["B8_flux_e"], self.syst_unc["Sur_prob_e"],
                                self.syst_unc["trigger_e"], self.syst_unc["reconstruction"]]
        HEP_CC_tot_syst_unc  = [self.syst_unc["HEP_flux_e"], self.syst_unc["Sur_prob_e"],
                                self.syst_unc["trigger_e"], self.syst_unc["reconstruction"]]
        HEP_ES_tot_syst_unc  = [self.syst_unc["HEP_flux_e"], self.syst_unc["Sur_prob_e"],
                               self.syst_unc["trigger_e"], self.syst_unc["reconstruction"]]
        Radon_tot_syst_unc   = [self.syst_unc["radon"], self.syst_unc["reconstruction"],
                                self.syst_unc["trigger_e"]]
        Neutron_tot_syst_unc = [self.syst_unc["neutron"], self.syst_unc["reconstruction"],
                                self.syst_unc["trigger_e"]]
        Ar42_tot_syst_unc    = [self.syst_unc["argon_upper"], self.syst_unc["reconstruction"],
                                self.syst_unc["trigger_e"]]

        O15_tot_syst_unc     = h.add_in_quadrature(O15_tot_syst_unc    )
        F17_tot_syst_unc     = h.add_in_quadrature(F17_tot_syst_unc    )
        B8_CC_tot_syst_unc   = h.add_in_quadrature(B8_CC_tot_syst_unc  )
        B8_ES_tot_syst_unc   = h.add_in_quadrature(B8_ES_tot_syst_unc  )
        HEP_CC_tot_syst_unc  = h.add_in_quadrature(HEP_CC_tot_syst_unc )
        HEP_ES_tot_syst_unc  = h.add_in_quadrature(HEP_ES_tot_syst_unc )
        Radon_tot_syst_unc   = h.add_in_quadrature(Radon_tot_syst_unc  )
        Neutron_tot_syst_unc = h.add_in_quadrature(Neutron_tot_syst_unc)
        Ar42_tot_syst_unc    = h.add_in_quadrature(Ar42_tot_syst_unc   )

        # print(O15_tot_syst_unc    )
        # print(F17_tot_syst_unc    )
        # print(B8_CC_tot_syst_unc  )
        # print(B8_ES_tot_syst_unc  )
        # print(HEP_CC_tot_syst_unc )
        # print(HEP_ES_tot_syst_unc )
        # print(Radon_tot_syst_unc  )
        # print(Neutron_tot_syst_unc)
        # print(Ar42_tot_syst_unc   )
        
        if self.config.metallicity == "high": module = IDh
        if self.config.metallicity == "low" : module = IDl
        
        for i in tqdm(range(self.ntests)):
            O15_test     = shuffle_slightly(module.O15    .map, O15_tot_syst_unc    , self.low, self.high)
            F17_test     = shuffle_slightly(module.F17    .map, F17_tot_syst_unc    , self.low, self.high)
            B8_CC_test   = shuffle_slightly(module.B8_CC  .map, B8_CC_tot_syst_unc  , self.low, self.high)
            B8_ES_test   = shuffle_slightly(module.B8_ES  .map, B8_ES_tot_syst_unc  , self.low, self.high)
            HEP_CC_test  = shuffle_slightly(module.HEP_CC .map, HEP_CC_tot_syst_unc , self.low, self.high)
            HEP_ES_test  = shuffle_slightly(module.HEP_ES .map, HEP_ES_tot_syst_unc , self.low, self.high)
            Radon_test   = shuffle_slightly(module.Radon  .map, Radon_tot_syst_unc  , self.low, self.high)
            Neutron_test = shuffle_slightly(module.Neutron.map, Neutron_tot_syst_unc, self.low, self.high)
            Ar42_test    = shuffle_slightly(module.Ar42   .map, Ar42_tot_syst_unc   , self.low, self.high)

            expected = {}
            observed = {}
            
            for i in O15_test: expected[i] = (B8_CC_test[i] + B8_ES_test[i] + HEP_CC_test[i] + HEP_ES_test[i] +
                                              Radon_test[i] + Neutron_test[i] + Ar42_test[i])

            for i in O15_test: observed[i] = (O15_test[i] + F17_test[i] + B8_CC_test[i] +
                                              B8_ES_test[i] + HEP_CC_test[i] + HEP_ES_test[i] +
                                              Radon_test[i] + Neutron_test[i] + Ar42_test[i])
            
            expected = np.array([expected[i] for i in expected])
            observed = np.array([observed[i] for i in observed])
            chi2_values.append(self.test._compute(expected, observed))
        return chi2_values
