from polpymer.core_funcs import Polymer, Monomer, Dish
from polpymer.data_funcs import plot_dish, plot_polymer, grow_polymer, \
     generate_N_polymers, expect_observ, error_observ
import numpy as np
import matplotlib.pyplot as plt



dish = Dish((10,10),(5,5))
dish.PERM(10, 2, 100)

plot_dish(dish, stems=False)

dish.polymer_correlation(bouqet=True)

plot_dish(dish, bouqet=True, stems=False)