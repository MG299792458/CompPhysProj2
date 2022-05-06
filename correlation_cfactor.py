from polpymer.core_funcs import Polymer, Monomer, Dish, correlation_metric
from polpymer.data_funcs import plot_dish, plot_polymer, grow_polymer, \
     generate_N_polymers, expect_observ, error_observ
import numpy as np
import matplotlib.pyplot as plt
from time import time

cfactors = [5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5]
metrics = []
times = []

for c in cfactors:
     for i in range(20):
          start = time()
          print('starting on: {} run: {}'.format(c,i))
          dish = Dish((10,10), (5,5))
          dish.PERM(10, c, 100)

          dish.polymer_correlation(bouqet=True)
          dish.correlation()

          metric = dish.corr_metric
          metrics.append(metric)

          dish = None
          end = time() - start
          times.append(end)
          print('finishing c: {} run: {}'.format(c, i))

cfactors = np.asarray(cfactors)
metrics = np.asarray(metrics)

np.save('Data/metrics', metrics)
cfactorsm = np.repeat(cfactors, 20)
print(cfactors)
print(metrics)

averages = []
deviations = []



plt.errorbar(cfactors, averages, yerr=deviations)
plt.scatter(cfactorsm, metrics, color='grey')
plt.xlabel(r'$c_{+} / c_{-}$')
plt.ylabel(r"Dish correlation $\tilde{\mathcal{C}}$")
plt.show()

