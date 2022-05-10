
from polpymer.core_funcs import Polymer, Monomer, Dish, correlation_metric
from polpymer.data_funcs import plot_dish, plot_polymer,\
      expect_observ, error_observ
from scipy.optimize import curve_fit as cv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from time import time

mpl.rcParams['figure.dpi'] = 180
mpl.rcParams['font.family'] = ["Baskerville"]
mpl.rcParams['font.size'] = 12

length = 70
cfac = 2

dish = Dish((10,10),(5,5))
start = time()
dish.PERM(100, cfac, length)
spent = time() - start


lengths = [polymer.chain_length for polymer in dish.polymers]
weights, observ = dish.weights, dish.end_to_end

e2e_weight_avg = expect_observ(observ, weights)
e2e_err = error_observ(observ, weights, 5)

def fun(xs, a, b):
    return a*xs**(2*b)

x = np.arange(length)
y = e2e_weight_avg

# idx = [ i for i in np.arange(40,len(e2e_err)) if np.isnan(e2e_err[i])or e2e_err[0] == 0]

# if len(idx) != 0:
#     cut_off = idx[0]
# else:
#     cut_off = -1

msk = [i for i in np.arange(len(x)) if e2e_err[i] != 0 and not np.isnan(e2e_err[i])]

copt, ccov = cv(fun, x[msk], y[msk], p0=[0.77,3/4], sigma=e2e_err[msk], absolute_sigma=False)

unc = np.sqrt(np.diag(ccov))

print(unc[1])

dish.correlation()
corr_metric = dish.corr_metric

exp_fit = fun(x[msk], *copt)
yerr = e2e_err

fig, ax = plt.subplots()
axx = ax.twinx()

print(len(x[msk]))

axx.plot(x, fun(x, *copt), label=r'fit $\alpha=$ '+'{:.2f}\t'.format(copt[0])+r"$\nu =$"+"{:.2f}".format(copt[1]), color='red', linestyle="--")
axx.errorbar(x[0:-1], e2e_weight_avg[0:-1], yerr=yerr[0:-1], label=r'exp', color='lightcoral', linestyle="None", marker="2")
axx.legend(frameon=False, loc=9)
axx.set_ylabel(r"$\langle r_e^2(L) \rangle$")
axx.set_xlabel(r"Polymer length $L$ [$a$]")


ax.hist(lengths, length, color='cornflowerblue')
ax.set_ylabel(r"counts [$-$]")

plt.title("End to End Distance PERM\n"+r"$\tilde{\mathcal{C}}=$"+"{:.2f}".format(corr_metric)+"\t"+r"$c_+=$"+"{}".format(cfac)+"\t"+r"$N=$"+"{}".format(200)+"\t"+r"$L=$"+"{}".format(length))
plt.gcf().set_size_inches(8,5)
ax.set_xlabel(r"polymer length L [$-$]")
plt.savefig("Figures/e2e_dist_P.pdf")
plt.show()

length = 50

dish = None
start = time()
dish = Dish((10,10),(5,5))
dish.find_N_polymer(800,length)
dish.analyse_polymers(length)
spent = time() - start


lengths = [polymer.chain_length for polymer in dish.polymers]
weights, observ = dish.weights, dish.end_to_end

e2e_weight_avg = expect_observ(observ, weights)
e2e_err = error_observ(observ, weights, 5)

def fun(xs, a, b):
    return a*xs**(2*b)

x = np.arange(length)
y = e2e_weight_avg


msk = [i for i in np.arange(len(x)) if e2e_err[i] != 0 and not np.isnan(e2e_err[i])]

copt, ccov = cv(fun, x[msk], y[msk], p0=[0.77,3/4], sigma=e2e_err[msk], absolute_sigma=False)

dish.correlation()
corr_metric = dish.corr_metric

exp_fit = fun(x[msk], *copt)
yerr = e2e_err

fig, ax = plt.subplots()
axx = ax.twinx()

print(len(x[msk]))

axx.plot(x, fun(x, *copt), label=r'fit $\alpha=$ '+'{:.2f}\t'.format(copt[0])+r"$\nu =$"+"{:.2f}".format(copt[1]), color='red', linestyle="--")
axx.errorbar(x[0:-1], e2e_weight_avg[0:-1], yerr=yerr[0:-1], label=r'exp', color='lightcoral', linestyle="None", marker="2")
axx.legend(frameon=False, loc=9)
axx.set_ylabel(r"$\langle r_e^2(L) \rangle$")
axx.set_xlabel(r"Polymer length $L$ [$a$]")


ax.hist(lengths, 50, color='cornflowerblue')
ax.set_ylabel(r"counts [$-$]")

plt.title("End to End Distance Rosenbluth\n"+r"$\tilde{\mathcal{C}}=$"+"{:.2f}".format(corr_metric)+"\t"+r"$N=$"+"{}".format(800)+"\t"+r"$L=$"+"{}".format(length))
plt.gcf().set_size_inches(8,5)
ax.set_xlabel(r"polymer length L [$-$]")
plt.savefig("Figures/e2e_dist_R.pdf")
plt.show()