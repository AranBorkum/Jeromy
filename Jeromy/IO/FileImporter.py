import random
import numpy as np
from collections import OrderedDict

_all = OrderedDict()


def get(*name):
    if len(name) == 0:
        return None
    if len(name) == 1:
        return _all[name[0]] if name[0] in _all else None

    return [_all[r] for r in name if name in _all]

def get_all():
    """Returns a list of all variables."""
    return _all.values()

class ImportFile():

    def __init__(self,
                 name = "",
                 path = "",
                 file = "",
                 xValues = [],
                 yValues = [],
                 map = {}):
        
        self._name = name
        self._path = path
        self._file_name = file
        self.map = map
        self.yValues = yValues
        self.xValues = xValues
        if path and path[-1] == "/":
            inFile = np.loadtxt(path + file, delimiter=",")
            self._file = inFile
            self.xValues = inFile[:, 0]
            self.yValues = inFile[:, 1]
            dataMap = {}
            for i, x in enumerate(inFile[:, 0]):
                dataMap[round(x, 2)] = inFile[:, 1][i]
            self.map = dataMap

                
        elif path:
            inFile = np.loadtxt(path + "/" + file, delimiter=",")
            self._file = inFile
            self.xValues = inFile[:, 0]
            self.yValues = inFile[:, 1]
            dataMap = {}
            for i, x in enumerate(inFile[:, 0]):
                dataMap[round(x, 2)] = inFile[:, 1][i]
            self.map = dataMap
        _all[name] = self

    def set_name(self, name):
        self._name = name

    # def combine_inputs(in1, in2):
    #     output = {}
    #     for i in in1.map:
    #         if i in in2.map:
    #             value = in1.map[i] + in2.map[i]
    #             output[i] = value

    #     return output

    def combine_inputs(self, list):

        xValues = [i for i in list[0].xValues]
        yValues = [0 for i in xValues]

        for item in list:
            for i, x in enumerate(xValues):
                if x in item.map:
                    yValues[i] += item.map[x]

        new_map = {}
        for i, x in enumerate(xValues):
            new_map[round(x, 1)] = yValues[i]
        self.xValues = xValues
        self.yValues = yValues
        self.map     = new_map
            
    def shuffle_slightly(self, error, low, high):
        new_file = ImportFile(name="{self.name}_shuffled")
        map = {}
        for i in self.map:
            if i >= low and i <= high:
                stat_error = np.sqrt(self.map[i])
                syst_error = self.map[i] * error
                map[i] = self.map[i] + (stat_error + syst_error)*random.randint(-1, 2)
                # seed = random.randint(-1, 2)
                # if seed > 0: map[i] = self.map[i]*1.05
                # if seed < 0: map[i] = self.map[i]*0.95
                # if seed == 0: map[i] = self.map[i]
        new_file.map = map
        new_file.xValues = [i for i in map]
        new_file.yValues = [map[i] for i in map]
        return new_file

        
    def scale(self, time_scale=1, background_scale=0):
        output = self
        new_map = {}
        yValues = []

        for i, y in enumerate(self.yValues):
            new_value = y * float(time_scale)
            if background_scale: new_value /= float(background_scale)            
            new_map[self.xValues[i]] = new_value
            yValues.append(new_value)
            
        output.map = new_map
        output.yValues = yValues
        return output


        
        
__all__ = [ImportFile]
