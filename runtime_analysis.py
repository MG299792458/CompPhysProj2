from polpymer.core_funcs import Polymer, Monomer, Dish, correlation_metric
from polpymer.data_funcs import plot_dish, plot_polymer, \
     expect_observ, error_observ
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from time import time
from scipy.optimize import curve_fit as cv

import cProfile

dish = Dish((10,10),(5,5))

cProfile.run("dish.PERM(100,2,70)", sort='tottime')