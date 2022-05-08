
from polpymer.core_funcs import Polymer, Monomer, Dish, correlation_metric
from polpymer.data_funcs import plot_dish, plot_polymer, grow_polymer, \
     generate_N_polymers, expect_observ, error_observ
from scipy.optimize import curve_fit as cv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 180
mpl.rcParams['font.family'] = ["Baskerville"]
mpl.rcParams['font.size'] = 12

length = 90
cfac = 2

dish = Dish((10,10),(5,5))
dish.PERM(1000, cfac, length)

lengths = [polymer.chain_length for polymer in dish.polymers]
weights, observ = dish.weights, dish.end_to_end

e2e_weight_avg = expect_observ(observ, weights)
e2e_err = error_observ(observ, weights, 5)

def fun(xs, a, b):
    return a*xs**(3/2)+b

x = np.arange(length)
y = e2e_weight_avg

idx = [ i for i in np.arange(len(e2e_weight_avg)) if np.isnan(e2e_weight_avg[i])]

if len(idx) != 0:
    cut_off = idx[0]
else:
    cut_off = -1

copt, ccov = cv(fun, x[0:cut_off], y[0:cut_off], sigma=e2e_err[0:cut_off], absolute_sigma=True)

dish.correlation()
corr_metric = dish.corr_metric

exp_fit = fun(x[0:cut_off], *copt)
yerr = e2e_err

fig, ax = plt.subplots()
axx = ax.twinx()

print(len(x[0:cut_off]))

axx.plot(x, fun(x, *copt), label=r'fit $\alpha=$ '+'{:.2f}'.format(copt[0]), color='red', linestyle="--")
axx.errorbar(x[0:-1], e2e_weight_avg[0:-1], yerr=yerr[0:-1], label=r'exp', color='lightcoral', linestyle="None", marker="2")
axx.legend(frameon=False, loc=9)
axx.set_ylabel(r"$\langle r_e^2(L) \rangle$")
axx.set_xlabel(r"Polymer length $L$ [$a$]")


ax.hist(lengths, 50, color='cornflowerblue')
ax.set_ylabel(r"counts [$-$]")

plt.title(r"End to End Distance"+"\t"+r"$\tilde{\mathcal{C}}=$"+"{:.2f}\t".format(corr_metric)+r"$c_+ / c_- =$"+"{}".format(cfac))
plt.gcf().set_size_inches(8,5)
plt.savefig("Figures/e2e_dist.pdf")
plt.show()

plt.plot(e2e_err)
plt.semilogy()
plt.show()