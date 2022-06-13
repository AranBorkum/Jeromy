import matplotlib.pyplot as plt

background_reduction = [1, 10, 20, 30, 50, 100, 200, 1000]
lowZ_standard      = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  -1.0, 1.02 ]
highZ_standard     = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  -1.0, 2.03 ]

lowZ_optimistic    = [0.02, 0.36, 0.81, 0.94, 1.93, 3.38,  5.82, 11.10]
highZ_optimistic   = [0.06, 0.82, 1.59, 2.41, 3.81, 6.76, 11.77, 20.39]

lowZ_better_boron  = [0.01, 0.38, 0.65, 1.18, 1.83, 3.57,  6.57, 15.53]
highZ_better_boron = [0.20, 0.69, 1.62, 2.33, 3.73, 6.70, 12.00, 22.36]

plt.figure()
plt.plot(background_reduction, highZ_standard    , "^", color="blue")
plt.plot(background_reduction, lowZ_standard     , "s", color="blue")

plt.plot(background_reduction, highZ_optimistic  , "^", color="orange")
plt.plot(background_reduction, lowZ_optimistic   , "s", color="orange")

plt.plot(background_reduction, highZ_better_boron, "^", color="green")
plt.plot(background_reduction, lowZ_better_boron , "s", color="green")

plt.plot([0, 1000], [5, 5], "k--")
plt.xscale("log")
plt.yscale("log")
plt.ylim([0, 25])
plt.show()

