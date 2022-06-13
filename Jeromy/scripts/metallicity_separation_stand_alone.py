#!/usr/bin/env python3

import Jeromy
import Jeromy.utilities.helpers as h

import numpy as np
import matplotlib.pyplot as plt

true_value_highZ = 250593.560
true_error_highZ = 65777.826

true_value_lowZ = 140384.420
true_error_lowZ = 35955.955

trials_highZ = {"standard"                                     : [250492.31, 654151152.92],
                "standard low background"                      : [250492.31, 920042.71   ],
                "optimistic low background"                    : [250492.31, 165071.07   ],
                "optimistic better boron"                      : [250492.31, 124505.26   ],
                "optimistic better boron tighter reco and xsec": [250492.31, 94204.09    ]}

trials_lowZ = {"standard"                                     : [140339.85, 593676748.23],
               "standard low background"                      : [140339.85, 793046.47   ],
               "optimistic low background"                    : [140339.85, 137455.27   ],
               "optimistic better boron"                      : [140339.85, 105522.06   ],
               "optimistic better boron tighter reco and xsec": [140339.85, 77854.38    ]}

if __name__ == "__main__":

    plt.figure()
    plt.subplot(1, 2, 1)
    x_values_high = np.linspace(true_value_highZ-(5*true_error_highZ),
                                true_value_highZ+(5*true_error_highZ),
                                1000)
    plt.plot(x_values_high, h.gaussian(x_values_high, true_value_highZ, true_error_highZ, 1))
    for i in trials_highZ:
        plt.plot(x_values_high, h.gaussian(x_values_high, trials_highZ[i][0], trials_highZ[i][1], 1), label=i)

    plt.legend(loc="best")

    plt.subplot(1, 2, 2)
    x_values_low = np.linspace(true_value_lowZ-(5*true_error_lowZ),
                                true_value_lowZ+(5*true_error_lowZ),
                                1000)
    plt.plot(x_values_low, h.gaussian(x_values_low, true_value_lowZ, true_error_lowZ, 1))
    for i in trials_lowZ:
        plt.plot(x_values_low, h.gaussian(x_values_low, trials_lowZ[i][0], trials_lowZ[i][1], 1), label=i)

    plt.legend(loc="best")
    plt.show()
    


    plt.figure()
    plt.subplot(1, 2, 1)
    plt.plot(x_values_high, h.gaussian(x_values_high, true_value_highZ, true_error_highZ, 1))
    for i in trials_lowZ:
        plt.plot(x_values_low, h.gaussian(x_values_low, trials_lowZ[i][0], trials_lowZ[i][1], 1), label=i)

    plt.subplot(1, 2, 2)
    plt.plot(x_values_low, h.gaussian(x_values_low, true_value_lowZ, true_error_lowZ, 1))
    for i in trials_highZ:
        plt.plot(x_values_high, h.gaussian(x_values_high, trials_highZ[i][0], trials_highZ[i][1], 1), label=i)


    plt.show()






