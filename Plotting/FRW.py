from polpymer.core_funcs import Polymer, Monomer, Dish, correlation_metric
from polpymer.data_funcs import plot_dish, plot_polymer,\
      expect_observ, error_observ
from scipy.optimize import curve_fit as cv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

dish = Dish((10,10),(5,5))
dish.FRW(10000,100)

end_to_end = dish.end_to_end
gyration = dish.gyration
end_to_end_avg = np.mean(end_to_end, 0)
gyration_avg = np.mean(gyration, 0)
errore2e = np.std(end_to_end, 0)/np.sqrt(10000)
errorgyr = np.std(gyration ,0)/np.sqrt(10000)

ls = np.arange(100)+1
plt.figure()
plt.errorbar(ls, end_to_end_avg, errore2e, label = 'End-to-end distance')
plt.xlabel('Polymer length')
plt.ylabel('End-to-end distance')
plt.savefig('.../Figures/FRWEnd2End.png')

plt.figure()
plt.errorbar(ls, gyration_avg, errorgyr, label = 'Radius of gyration')
plt.xlabel('Polymer length')
plt.ylabel('Radius of gyration')
plt.savefig('.../Figures/FRWRadiusofGyration.png')