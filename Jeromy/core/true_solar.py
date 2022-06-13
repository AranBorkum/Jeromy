import yaml
import math
import Jeromy
import Jeromy.utilities.helpers as h
import Jeromy.IO.UncertaintyImporter as UncertaintyImporter
import Jeromy.config.configuration_loader   as configuration_loader
import numpy as np

from scipy.integrate import quad 
import matplotlib.pyplot as plt

def gaussian(x, m, s, a):
    term1 = a / (s * np.sqrt(2*np.pi))
    term2 = (x-m)**2
    term3 = 2 * (s**2)
    return term1 * np.exp(-term2 / term3)

def orderOfMagnitude(number):
    return math.floor(math.log(number, 10))

with open(Jeromy.__jeromy_data__+"/solar_rates.yml", "r") as stream:
    solar_rates = yaml.safe_load(stream)

class _spectrum():

    def __init__(self, name, mean, error):
        self.name = name
        self.mean = mean
        self.error = error

        order = orderOfMagnitude(self.error)
        self.constant = 10 ** order
        
        self.parameters = [self.mean, self.error, 1]
        self.function = gaussian

    def p_value(self, test_stat):

        if test_stat > self.mean:
            I, e = quad(gaussian, test_stat/self.constant, np.inf,
                        args=(self.mean/self.constant, self.error/self.constant, 1))
        else:
            I, e = quad(gaussian, -np.inf, test_stat/self.constant,
                        args=(self.mean/self.constant, self.error/self.constant, 1))

        return 2*I

    def plot(self):
        plt.figure()
        xValues = np.linspace(self.mean - (self.error*5), self.mean + (self.error*5), 1000)
        yValues = [self.function(i, *self.parameters) for i in xValues]
        plt.plot(xValues, yValues)
        
        plt.plot([self.mean, self.mean], [0, self.function(self.mean, *self.parameters)], "k--")
        plt.plot([self.mean+self.error, self.mean+self.error],
                 [0, self.function(self.mean+self.error, *self.parameters)], "k--")
        plt.plot([self.mean-self.error, self.mean-self.error],
                 [0, self.function(self.mean-self.error, *self.parameters)], "k--")
        plt.plot([self.mean+self.error, self.mean-self.error],
                 [self.function(self.mean-self.error, *self.parameters),
                  self.function(self.mean-self.error, *self.parameters)], "k--")

        plt.plot([self.mean+(self.error*2), self.mean+(self.error*2)],
                 [0, self.function(self.mean+(self.error*2), *self.parameters)], "k--")
        plt.plot([self.mean-(self.error*2), self.mean-(self.error*2)],
                 [0, self.function(self.mean-(self.error*2), *self.parameters)], "k--")
        plt.plot([self.mean+(self.error*2), self.mean-(self.error*2)],
                 [self.function(self.mean-(self.error*2), *self.parameters),
                  self.function(self.mean-(self.error*2), *self.parameters)], "k--")

        plt.plot([self.mean+(self.error*3), self.mean+(self.error*3)],
                 [0, self.function(self.mean+(self.error*3), *self.parameters)], "k--")
        plt.plot([self.mean-(self.error*3), self.mean-(self.error*3)],
                 [0, self.function(self.mean-(self.error*3), *self.parameters)], "k--")
        plt.plot([self.mean+(self.error*3), self.mean-(self.error*3)],
                 [self.function(self.mean-(self.error*3), *self.parameters),
                  self.function(self.mean-(self.error*3), *self.parameters)], "k--")

        plt.xlabel("Rate")
        plt.ylabel("Probability")
        plt.show()

        
class SolarTrue():

    def __init__(self, config):
        self.config = config

        config_name = self.config.name
        if config_name[:3] == "low":
            conf_high = configuration_loader.ConfigurationLoader("high" + self.config.name[3:])
            conf_low  = configuration_loader.ConfigurationLoader("low"  + self.config.name[3:])
            syst_high = conf_high.systematic_uncertainties
            syst_low  = conf_low .systematic_uncertainties
        else:
            conf_high = configuration_loader.ConfigurationLoader("high" + self.config.name[4:])
            conf_low  = configuration_loader.ConfigurationLoader("low"  + self.config.name[4:])
            syst_high = conf_high.systematic_uncertainties
            syst_low  = conf_low .systematic_uncertainties

        
        O15_rate_high = solar_rates["high"]["O15"]/10*float(self.config.exposure)
        F17_rate_high = solar_rates["high"]["F17"]/10*float(self.config.exposure)
        CNO_rate_high = O15_rate_high + F17_rate_high

        O15_rate_low = solar_rates["low"]["O15"]/10*float(self.config.exposure)
        F17_rate_low = solar_rates["low"]["F17"]/10*float(self.config.exposure)
        CNO_rate_low = O15_rate_low + F17_rate_low

        O15_error_high = [syst_high["O15_flux_e"],
                          syst_high["Sur_prob_e"],
                          syst_high["trigger_e"]]
        F17_error_high = [syst_high["F17_flux_e"],
                          syst_high["Sur_prob_e"],
                          syst_high["trigger_e"]]
        O15_error_high = h.add_in_quadrature(O15_error_high)*O15_rate_high
        F17_error_high = h.add_in_quadrature(F17_error_high)*F17_rate_high
        CNO_error_high = h.add_in_quadrature([O15_error_high, F17_error_high])
        
        O15_error_low = [syst_low["O15_flux_e"],
                         syst_low["Sur_prob_e"],
                         syst_low["trigger_e"]]
        F17_error_low = [syst_low["F17_flux_e"],
                         syst_low["Sur_prob_e"],
                         syst_low["trigger_e"]]
        O15_error_low = h.add_in_quadrature(O15_error_low)*O15_rate_low
        F17_error_low = h.add_in_quadrature(F17_error_low)*F17_rate_low
        CNO_error_low = h.add_in_quadrature([O15_error_low, F17_error_low])
        
        self.O15_high = _spectrum("O15_high", O15_rate_high, O15_error_high)
        self.F17_high = _spectrum("F17_high", F17_rate_high, F17_error_high)
        self.CNO_high = _spectrum("CNO_high", CNO_rate_high, CNO_error_high)

        self.O15_low = _spectrum("O15_low", O15_rate_low, O15_error_low)
        self.F17_low = _spectrum("F17_low", F17_rate_low, F17_error_low)
        self.CNO_low = _spectrum("CNO_low", CNO_rate_low, CNO_error_low)



        
