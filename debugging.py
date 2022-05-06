from polpymer.core_funcs import Polymer, Monomer, Dish, correlation_metric
from polpymer.data_funcs import plot_dish, plot_polymer, grow_polymer, \
     generate_N_polymers, expect_observ, error_observ
import numpy as np
import matplotlib.pyplot as plt



dish = Dish((10,10),(5,5))
dish1 = Dish((10,10),(5,5))
dish.PERM(10, 10, 100)
dish1.PERM(10, 5, 100)

dish.polymer_correlation(bouqet=True)
dish1.polymer_correlation(bouqet=True)

plot_dish(dish, bouqet=True, stems=False)
plot_dish(dish1, bouqet=True)

dish.correlation()
dish1.correlation()

corr_matrix = dish.correlation_matrix
corr_metric = dish.corr_metric

corr_matrix1 = dish1.correlation_matrix
corr_metric1 = dish1.corr_metric



plt.imshow(corr_matrix, origin='upper')
plt.colorbar()
plt.xlabel(r"polymer $j$")
plt.ylabel(r"polymer $i$")
plt.title(r"correlation $r_{x}^{ij} \cdot r_{y}^{ij}$")
plt.scatter(18, 23, marker='.', s=20, color='red')
plt.show()


plt.imshow(corr_matrix1, origin='upper')
plt.colorbar()
plt.xlabel(r"polymer $j$")
plt.ylabel(r"polymer $i$")
plt.title(r"correlation $r_{x}^{ij} \cdot r_{y}^{ij}$")
plt.scatter(18, 23, marker='.', s=20, color='red')
plt.show()

print(corr_metric, corr_metric1)