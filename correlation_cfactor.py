from polpymer.core_funcs import Polymer, Monomer, Dish, correlation_metric
from polpymer.data_funcs import plot_dish, plot_polymer, grow_polymer, \
     generate_N_polymers, expect_observ, error_observ
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from time import time

mpl.rcParams['figure.dpi'] = 180
mpl.rcParams['font.family'] = ["Baskerville"]
mpl.rcParams['font.size'] = 12

realisations = 100
cfactors = [5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12]
metrics = []
times = []
baselines = []
baseline_times = []

for i in range(realisations):
     dish = Dish((10,10),(5,5))
     start = time()
     dish.find_N_polymer(10, 100)
     dish.polymer_correlation(bouqet=True)
     dish.correlation()
     end = time() - start
     baseline_times.append(end)
     baseline = dish.corr_metric
     baselines.append(baseline)
     dish = None

baseline = np.average(baselines)
baseline_std = np.std(baselines)

for c in cfactors:
     for i in range(realisations):
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
cfactorsm = np.repeat(cfactors, realisations)

averages = np.nanmean(metrics.reshape((len(cfactors), realisations)), axis=1)
deviations = np.nanstd(metrics.reshape((len(cfactors), realisations)), axis=1)

plt.axhline(baseline, color='blue', linestyle='--')
plt.fill_between(cfactors, baseline+baseline_std, baseline-baseline_std, color='powderblue')
plt.scatter(cfactorsm, metrics, color='grey', s=10)
plt.errorbar(cfactors, averages, yerr=deviations, color='red', linestyle="None", marker='s')
plt.xlabel(r'$c_{+} / c_{-}$ [$-$]')
plt.ylabel(r"$\tilde{\mathcal{C}}$ [-]")
plt.title(r'Correlation of Polymers')
plt.gcf().set_size_inches(8,5)
plt.savefig("Figures/correlation_per_cfac.pdf")
plt.show()

baseline_time = np.average(baseline_times)
baseline_time_std = np.std(baseline_times)

times = np.asarray(times).reshape((len(cfactors), realisations))
avg_time = np.average(times, axis=1)
std_time = np.std(times, axis=1)

plt.axhline(baseline_time, color='blue', linestyle='--')
plt.fill_between(cfactors, baseline_time+baseline_time_std, baseline_time-baseline_time_std, color='powderblue')
plt.scatter(cfactorsm, times, color='grey', s=10)
plt.errorbar(cfactors, avg_time, yerr=std_time, color='red', linestyle='None', marker='s')
plt.xlabel(r"$c_{+}/c_{-}$ [$-$]")
plt.ylabel(r'$t_{elapsed}$ [$s$]')
plt.title(r'Time spent for Polymer Creation')
plt.gcf().set_size_inches(8,5)
plt.savefig("Figures/time_spent_cfac.pdf")
plt.show()
