from polpymer.core_funcs import Polymer, Monomer, Dish, correlation_metric
from polpymer.data_funcs import plot_dish, plot_polymer, \
     expect_observ, error_observ
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from time import time
from scipy.optimize import curve_fit as cv

mpl.rcParams['figure.dpi'] = 180
mpl.rcParams['font.family'] = ["Baskerville"]
mpl.rcParams['font.size'] = 12

realisations = 100
length = 50
cfactors = [2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8]
metrics = []

times = []
baselines = []
baseline_times = []

baseline_dist = []
dists = []

bins = []

amnts_base = []
amnts = []


for i in range(10):
     dish = Dish((10,10),(5,5))
     start = time()
     dish.find_N_polymer(800, length)
     end = time() - start
     dish.polymer_correlation(bouqet=True)
     dish.correlation()
     lengths = [polymer.chain_length for polymer in dish.polymers]
     amnts_base.append(len(lengths))
     dist, bins, bar = plt.hist(lengths, length)
     baseline_dist.append(dist)
     baseline_times.append(end)
     plt.close()
     baseline = dish.corr_metric
     print(str(i)+"\t Rosenbluth")
     baselines.append(baseline)
     dish = None

baseline = np.average(baselines)
baseline_std = np.std(baselines)

for c in cfactors:
     for i in range(realisations):
          print('starting on: {} run: {}'.format(c,i))
          dish = Dish((10,10), (5,5))
          start = time()
          dish.PERM(10, c, length)
          end = time() - start

          dish.polymer_correlation(bouqet=True)
          dish.correlation()

          metric = dish.corr_metric

          lengths = [polymer.chain_length for polymer in dish.polymers]
          amnts.append(len(lengths))
          dist, bins, bar = plt.hist(lengths, length)
          dists.append(dist)
          metrics.append(metric)
          plt.close()
          dish = None
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
plt.xlabel(r'$c_{+}$ [$-$]')
plt.ylabel(r"$\tilde{\mathcal{C}}$ [-]")
plt.title(r'Correlation of Polymers')
plt.subplots_adjust(left=0.09, right=0.9, bottom=0.09)
plt.gcf().set_size_inches(8,5)
plt.savefig("Figures/correlation_per_cfac.pdf")
plt.show()

baseline_time = np.average(baseline_times)
baseline_time_std = np.std(baseline_times)

times_n = np.asarray(times).reshape((len(cfactors), realisations))
avg_time = np.average(times_n, axis=1)
std_time = np.std(times_n, axis=1)

plt.axhline(baseline_time, color='blue', linestyle='--')
plt.fill_between(cfactors, baseline_time+baseline_time_std, baseline_time-baseline_time_std, color='powderblue')
plt.scatter(cfactorsm, times, color='grey', s=10)
plt.errorbar(cfactors, avg_time, yerr=std_time, color='red', linestyle='None', marker='s')
plt.xlabel(r"$c_{+}$ [$-$]")
plt.ylabel(r'$t_{elapsed}$ [$s$]')
plt.title(r'Time spent for Polymer Creation')
plt.subplots_adjust(left=0.09, right=0.9, bottom=0.09)
plt.gcf().set_size_inches(8,5)
plt.savefig("Figures/time_spent_cfac.pdf")
plt.show()


corr_rel_err = (metrics - baseline)
speedup = -(times_n - baseline_time) / baseline_time

avg_rel_err = (averages - baseline)
avg_rel_err_err = (deviations ) / baseline

avg_speedup = -(avg_time - baseline_time) / baseline_time
avg_time_err = (std_time) /baseline_time

for i in range(len(cfactors)):
     start = i*realisations
     end = (i+1)*realisations
     col = 1-(i+1)/len(cfactors)
     plt.scatter(speedup[i], corr_rel_err[start:end], color=str(col), s=10, label=str(cfactors[i]))
plt.errorbar(avg_speedup, avg_rel_err, yerr=avg_rel_err_err, xerr=avg_time_err, color='red', marker="s", linestyle="None")
plt.ylabel(r"$\Delta \tilde{\mathcal{C}}$ [$-$]")
plt.xlabel(r'$speedup$ [-]')
plt.title(r'Correlation increase, speedup tradeof')
plt.semilogx()
plt.subplots_adjust(left=0.09, right=0.9, bottom=0.09)
plt.legend(frameon=False)
plt.gcf().set_size_inches(8,5)
plt.savefig("Figures/trade_off.pdf")
plt.show()

baseline_dist = np.asarray(baseline_dist)
bas_dist_avg = np.average(baseline_dist, axis=0)
bas_dist_std = np.std(baseline_dist, axis=0)

dists = np.asarray(dists)
dists = dists.reshape((realisations,len(cfactors),length))

dists_avg = np.mean(dists, axis=0)
dists_std = np.std(dists, axis=0)

bins = np.arange(length)+1

# plt.scatter(bins, bas_dist_avg, color='blue', label="Ros. baseline.", marker="2")
# plt.fill_between(bins, bas_dist_avg+bas_dist_std, bas_dist_avg-bas_dist_std, color='powderblue')
# for i in range(dists_avg.shape[0]):
#      plt.errorbar(bins, dists_avg[i], yerr=dists_std[i], marker="2", label=str(cfactors[i]), linestyle="None")
# plt.xlabel(r"Polymer Length [$-$]")
# plt.ylabel(r"Count [$-$]")
# plt.subplots_adjust(left=0.09, right=0.9, bottom=0.09)
# plt.legend(frameon=False)
# plt.gcf().set_size_inches(8,5)
# plt.savefig("Figures/trade_off.pdf")
# plt.show()

amnts_base = np.asarray(amnts_base)
amnts_base_avg = np.mean(amnts_base)
amnts_base_std = np.std(amnts_base)

amnts = np.asarray(amnts)
amnts = amnts.reshape((len(cfactors),realisations))
amnts_avg = np.mean(amnts, axis=1)
amnts_std = np.std(amnts, axis=1)


plt.axhline(amnts_base_avg, color='blue', linestyle='--', marker="s")
plt.fill_between(cfactors, amnts_base_avg+amnts_base_std, amnts_base_avg-amnts_base_std, color='powderblue')
plt.scatter(cfactorsm, amnts, color='grey', s=10)
plt.errorbar(cfactors, amnts_avg, yerr=amnts_std, color='red', linestyle="None", marker='s')
bottom, top = plt.ylim()
plt.ylim([0,top])
plt.xlabel(r'$c_{+}$ [$-$]')
plt.ylabel(r"count [-]")
plt.title(r'Number of polymers for given $c_+$')
plt.semilogy()
plt.subplots_adjust(left=0.09, right=0.9, bottom=0.09)
plt.gcf().set_size_inches(8,5)
plt.savefig("Figures/amnt_per_cfac.pdf")
plt.show()