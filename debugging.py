from polpymer.core_funcs import Polymer, Monomer, Dish
from polpymer.data_funcs import plot_polymer, grow_polymer, \
     generate_N_polymers, expect_observ, error_observ
import numpy as np
import matplotlib.pyplot as plt



dish = Dish((10,10),(5,5))
dish.PERM(10, 2, 100)
print(len(dish.polymers))