import random
import numpy as np

class ToyMonteCarlo():
    def __init__(self, name=""):
        self._name = name

    def run_toy_mc(input_map, error_syst, low, high):
        output = {}
        for i in input_map:
            if i >= low and i <= high:
                stat_error = np.sqrt(input_map[i])
                syst_error = input_map[i] * error_syst
                output[i] = input_map[i] + np.sqrt(stat_error**2 + syst_error**2) * random.randint(-1, 2)
        return output

