# Simulating Polymers as self avoiding random walks
The contents of this repository allow for the simulation of polymers as self avoiding random walks, to this end Dish, Polomer and Monomer classes are avaiable to control the growth of polymers to a high degree as well as analyze the ensembles using methods of these classes

## Table of contents
[[_TOC_]]

## Installation
The module can be installed using pip or can be imported when the working file is placed in the same directory.

```python
from polpymer.data_funcs import *
from polpymer.core_funcs import *
```
to access the necessary functions and classes.

## Usage

```python
from polpymer.core_funcs import Polymer, Monomer, Dish, correlation_metric
from polpymer.data_funcs import plot_dish, plot_polymer,\
      expect_observ, error_observ
import matplotlib.pyplot as plt


length = 70 # Length to grow to
cfac = 2    # C_+ factor

dish = Dish((10,10),(5,5))      # Initialise a Petri dish
dish.PERM(100, cfac, length)    # Use PERM to generate polymers

weights, observ = dish.weights, dish.end_to_end # Assignment of the
                                                # automatically calculated
                                                # weights and observables

e2e_weight_avg = expect_observ(observ, weights) # Calculation of expectation
e2e_err = error_observ(observ, weights, 5)      # Calculation of error using
                                                # bootstrapping

dish.correlation()                  # Calculation of the correlation metric
corr_metric = dish.corr_metric      # Requesting value

lengths = [pol.chain_length for pol in dish.polymers]   # Easily get the lengths
                                                        # of the polymers using
                                                        # built-in methods

# Plotting of the results
plt.plot(lengths, e2e_weight_avg)
plt.show()
```

## Contributors

@sangersjeroen
@agefrancke2*
@mwglorie*

*On gitlab.kwant-project.org
