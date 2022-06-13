import matplotlib.pyplot as plt
import numpy as np

def gaussian(x, m, s, a):
    term1 = a / (s * np.sqrt(2*np.pi))
    term2 = (x-m)**2
    term3 = 2 * (s**2)
    return term1 * np.exp(-term2 / term3)


    

mean = 0.47465565798554643
sigma = 0.0015058972209989197

optimisitic_unc_br10_high  = 0.47083416378469595
better_boron_unc_br10_high = 0.4706837365845946

x_values = np.linspace(mean - 5*sigma, mean + 5*sigma, 1000)
y_values = [gaussian(i, mean, sigma, 1) for i in x_values]

def plot_line(value):
    plt.plot([value, value], [0, 1.1*max(y_values)], "--")


plt.figure()
plt.plot(x_values, y_values)
plot_line(optimisitic_unc_br10_high)
plot_line(better_boron_unc_br10_high)
plt.show()
