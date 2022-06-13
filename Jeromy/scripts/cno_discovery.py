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
    systematics = config.systematic_uncertaintiesx
    
    # Read in the data for the given metallicity
    if config.metallicity = "high": dataIO = DataHZ
    if config.metallicity = "low" : dataIO = DataLZ

    





    
