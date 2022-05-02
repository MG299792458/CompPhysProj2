from polpymer.core_funcs import Polymer, Monomer, Dish
from polpymer.data_funcs import plot_polymer, grow_polymer, \
     generate_N_polymers, expect_observ, error_observ
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 120

dish = Dish((10,10),(5,5))
dish.PERM(100, 10, 4)


lengths = []
print(dish.polymers)

for polymer in dish.polymers:
    # plot_polymer(polymer)
    lengths.append(polymer.chain_length)


plt.hist(lengths)
plt.show()
