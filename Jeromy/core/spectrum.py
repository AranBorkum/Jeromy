import numpy as np

class Spectrum()

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


    # def combine_inputs(in1, in2):
    #     output = {}
    #     for i in in1.map:
    #         if i in in2.map:
    #             value = in1.map[i] + in2.map[i]
    #             output[i] = value

    #     return output

    def combine_inputs(list):
        output = {}
        for i in list[0].map:
            value = list[0].map[i]
            for item in list[1:]:
                if i in item.map:
                    value += item.map[i]
            output[i] = value
        return output
