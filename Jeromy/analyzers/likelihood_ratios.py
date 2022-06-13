import math
import yaml
import Jeromy
import Jeromy.core.true_solar as true_solar
import Jeromy.analyzers.counting_experiment as counting_experiment
import Jeromy.statistics.test_statistic as test_statistic
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from math import isnan
from lmfit.models import SkewedGaussianModel, GaussianModel

def gaussian(x, m, s, a):
    term1 = a / (s * np.sqrt(2*np.pi))
    term2 = (x-m)**2
    term3 = 2 * (s**2)
    return term1 * np.exp(-term2 / term3)

with open(Jeromy.__jeromy_data__+"/solar_rates.yml", "r") as stream:
    solar_rates = yaml.safe_load(stream)

def standard_deviation(mean, values):
    std_dev = 0
    for i in values: std_dev += (i-mean)**2
    std_dev /= len(values)
    return np.sqrt(std_dev)
    
class LikelihoodTest():

    def __init__(self, observed, test, config, anti):
        self._observed  = observed
        self.test       = test
        self.config     = config
        self.anti       = anti
        self.observed_mean  = self._observed.cno_rate
        self.observed_sigma = self._observed.cno_uncertainty

        solar_true = true_solar.SolarTrue(self.config)
        if config.metallicity == "high" and not anti:
            self._predicted = solar_true.CNO_high
        if config.metallicity == "high" and anti:
            self._predicted = solar_true.CNO_low
        if config.metallicity == "low" and not anti:
            self._predicted = solar_true.CNO_low
        if config.metallicity == "low" and anti:
            self._predicted = solar_true.CNO_high

        self.bins      = self.pre_process()[0]
        self.observed  = self.pre_process()[1]
        self.predicted = self.pre_process()[2]
        self.chi2      = self.test_stat(self.observed, self.predicted)
        self.true_bins = self.true_predicted_spectrum()[0]
        self.true_pred = self.true_predicted_spectrum()[1]

    def test_stat(self, obs, pre):
        test_statistic = 0
        for i in range(len(self.observed)):
            if pre[i]:
                test_statistic += pre[i] - obs[i] + obs[i] * np.log(obs[i] / pre[i])

        return test_statistic * 2
    
        
    def pre_process(self):
        observed  = []
        predicted = []

        bins = np.linspace(self._observed.cno_rate-(5*self._observed.cno_uncertainty),
                           self._observed.cno_rate+(5*self._observed.cno_uncertainty),
                           50)
        for i in bins:
            observed .append(gaussian(i, self._observed.cno_rate, self._observed.cno_uncertainty, 1))
            predicted.append(gaussian(i, self._predicted.mean, self._predicted.error, 1))

        return bins, observed, predicted

    def true_predicted_spectrum(self):
        bins = np.linspace(self._predicted.mean - (5*self._predicted.error),
                           self._predicted.mean + (5*self._predicted.error),
                           50)
        return bins, [gaussian(i, self._predicted.mean, self._predicted.error, 1) for i in bins]

    
        
class LikelihoodTestEnsemble():

    def __init__(self, name, roi, config, ntrials=1000, anti=False):
        self.name    = name
        self.roi     = roi
        self.config  = config
        self.ntrials = ntrials
        self.anti    = anti

        chi2_values = []
        for i in tqdm(range(self.ntrials)):
            trial = counting_experiment.CountingExperiment(self.name, self.roi, self.config, True)
            likelihood = LikelihoodTest(trial, None, config, self.anti)
            if not math.isinf(likelihood.chi2) and not math.isnan(likelihood.chi2):
                chi2_values.append(likelihood.chi2)
        self.chi2_values = chi2_values

    def plot_test_statistic(self, bins=50, title=""):
        plt.figure()
        plt.hist(self.chi2_values, bins=bins)
        plt.xlabel(r"Test statistic, $q_{0}$", fontsize=16)
        plt.title(title)
        plt.show()
        
    def make_test_statistic_histogram(self, bins=50, fit=False, model=None):
        hist = plt.hist(self.chi2_values, bins=bins)
        if fit and model:
            if model == "GaussianModel":
                model = GaussianModel()
            if model == "SkewedGaussianModel":
                model = SkewedGaussianModel()
            # find the mean of the histogram
            mean = 0
            bin_content = 0
            for i, x in enumerate(hist[1][:-1]):
                if hist[0][i] > bin_content:
                    bin_content = hist[0][i]
                    mean = x
            sigma = standard_deviation(mean, self.chi2_values)
            params = model.make_params(amplitude=bin_content, center=mean, sigma=sigma)
            result = model.fit(hist[0], params, x=hist[1][:-1])
            my_fit = plt.plot(hist[1][:-1], result.best_fit)

            
            
        return hist, my_fit
