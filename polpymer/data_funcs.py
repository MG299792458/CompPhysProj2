""" Polpymer Core Functions

This script contains all the functions used to simulate the creation of random
polymers by using a Monte-Carlo method.

Results of the simulation can be analysed using the functions in data_funcs.py
"""


# Module imports
from inspect import trace
import matplotlib.pyplot as plt
from polpymer.core_funcs import Polymer, Monomer, Dish
from random import choice
import numpy as np



# Global variables


# Testing for file import
if __name__ == "__main__":
    print("Please import functions from this file by using \
     'from polpymer import *', instead of running this file directly")


# Functions
def plot_polymer(polymer: object) -> None:

    x_ = np.asarray([])
    y_ = np.asarray([])

    xlines = np.asarray([])
    xlines_posx = np.asarray([])
    xlines_posy = np.asarray([])
    ylines = np.asarray([])
    ylines_posy = np.asarray([])
    ylines_posx = np.asarray([])

    (xmax, ymax) = polymer.dimensions
    cnt = 0

    for monomer in polymer:
        start = monomer.location
        end = monomer.end_location

        x_ = np.append(x_, start[0])
        y_ = np.append(y_, start[1])

        plt.text(x_[-1]+0.2, y_[-1]+0.2, str(cnt) )
        cnt += 1
        ang = monomer.angle
        if ang == 0 or ang == 2:
            if ang == 0:
                xlines_posx = np.append(xlines_posx, start[0])
            else:
                xlines_posx = np.append(xlines_posx, end[0])
            xlines = np.append(xlines, 1)
            xlines_posy = np.append(xlines_posy, start[1])
        if ang == 1 or ang == 3:
            if ang == 1:
                ylines_posy = np.append(ylines_posy, start[1])
            else:
                ylines_posy = np.append(ylines_posy, end[1])
            ylines = np.append(ylines, 1)
            ylines_posx = np.append(ylines_posx, start[0])

    plt.xlim([0,xmax])
    plt.ylim([0,ymax])

    plt.scatter(x_, y_, linestyle='None', marker='o', s=40)
    plt.vlines(ylines_posx, ylines_posy, ylines_posy+1, color='black')
    plt.hlines(xlines_posy, xlines_posx, xlines_posx+1, color='black')
    plt.show()


def plot_dish(dish: object, bouqet: bool=False, stems: bool=False) -> None:

    if bouqet and dish.bouqet is None:
        dish.polymer_correlation(bouqet=True)
    elif bouqet and dish.bouqet is not None:
        polymers = dish.bouqet
    else:
        polymers = dish.polymers

    if stems:
        x_ = np.asarray([])
        y_ = np.asarray([])

        xlines = np.asarray([])
        xlines_posx = np.asarray([])
        xlines_posy = np.asarray([])
        ylines = np.asarray([])
        ylines_posy = np.asarray([])
        ylines_posx = np.asarray([])

        for polymer in polymers:
            for monomer in polymer:
                start = monomer.location
                end = monomer.end_location

                x_ = np.append(x_, start[0])
                y_ = np.append(y_, start[1])

                ang = monomer.angle
                if ang == 0 or ang == 2:
                    if ang == 0:
                        xlines_posx = np.append(xlines_posx, start[0])
                    else:
                        xlines_posx = np.append(xlines_posx, end[0])
                    xlines = np.append(xlines, 1)
                    xlines_posy = np.append(xlines_posy, start[1])
                if ang == 1 or ang == 3:
                    if ang == 1:
                        ylines_posy = np.append(ylines_posy, start[1])
                    else:
                        ylines_posy = np.append(ylines_posy, end[1])
                    ylines = np.append(ylines, 1)
                    ylines_posx = np.append(ylines_posx, start[0])

        plt.vlines(ylines_posx, ylines_posy, ylines_posy+1, color='black')
        plt.hlines(xlines_posy, xlines_posx, xlines_posx+1, color='black')

    all_coords: list[tuple[int,int], ...] = []

    # Add all crossed coordinates to all_coords list
    for polymer in polymers:
        for monomer in polymer:
            coord_end = monomer.end_location
            all_coords.append(coord_end)

    # Parse the all_coords list for number of occurences
    occur: dict[tuple[int,int], int] = {}
    for coords in all_coords:
        if coords not in occur:
            occur[coords] = 1
        elif coords in occur:
            occur[coords] += 1

    # Process data in plottable sets
    xs = np.asarray([])
    ys = np.asarray([])
    ss = np.asarray([])
    for key in occur:
        xs = np.append(xs, key[0])
        ys = np.append(ys, key[1])
        ss = np.append(ss, occur[key])

    plt.scatter(xs, ys, ss)
    plt.show()


def expect_observ(observ, w):
    """function to calculate the expectation value for the observable

    Parameter
    ---------
    observ : nd.array
        the observable for which the expectation value has to be determined
    w : nd.array
        the weight of the polymer

    Return
    ------
    expect : nd.array
        expectation value of the observable
    """

    step_1 = w * observ
    step_2a = np.sum(step_1, axis=0)
    step_2b = np.sum(w, axis=0)

    expect = step_2a / step_2b

    return expect


def error_observ(observ, w, n):
    """function to calculate the error of the expectation value of the observable

    Parameter
    ---------
    observ : nd.array

    w : nd.array

    n : int

    Return
    ------
    error : nd.array
    """

    N, L = np.shape(observ)

    expect = np.zeros((n,L))

    for i in range(n):
        rand_sample = np.int_(np.round(np.random.rand(N,) * (N-0.5), decimals=0))

        rand_observ = observ[rand_sample, :]
        rand_w = w[rand_sample, :]

        expect_i = expect_observ(rand_observ, rand_w)
        expect[i,:] = expect_i

    error = np.sqrt(1/n * np.sum(expect**2, axis=0) - (1/n * np.sum(expect, axis=0))**2)

    return error



def scale_law(x, a, b):
    """Random walk scaling law.
    For the free random walk the value of b=1/2 and for the self avoiding random walk
    the value for b=3/4.
    
    Parameter
    ---------
    x : nd.array
        Length of the polymer
    a : float
        Amplitude, for the end-to-end distance of a 2D self avoiding random walk 0.771
    b : float
        Critical exponent
    
    Return
    ------
    y : nd.array
        Expectation value of the observable, like end-to-end distance and radius of gyration,
        at length x.
    """
    
    y = a * x**(2 * b)
    return y



def fit_observ_error(observ, error_observ, scale_val=np.array([0.771, 3/4])):
    """Automatically fits the simulated observable using the error of the observable and then 
    plot the simulated observable (with errorbars) against the fit and the expected behaviour
    of the observable based on the scaling law.
    
    Parameter
    ---------
    observ : nd.array
        The weighted averages of the observable for different polymer lengths
    error_observ : nd.array
        The error of the weighted averages of th observable for different polymer lengths
    scale_val: nd.array
        The expected values for the scaling factor according to literature. Must be an vector
        with two ellements ([a, b]), with the element relating to the scaling law as:
        <observ(L)> = a * L**(2 * b)
    
    Return
    ------
    length : nd.array
        The length of the polymer from 1 to L
    copt : nd.array
        The fitted values for scale_val with copt[0]=a and copt[1]=b
    ccov : nd.array
        The covariance of the fitted values
    expect_fit : nd.array
        The expected weighted average of the observable based on scale_val
    fit_plot_val: nd.array
        The fit of the weighted average of the observable based on copt
    """
    
    length = np.arange(1, len(observ)+1, 1)
    
    copt, ccov = cv(scale_law, length[1:], observ[1:], p0=scale_val,
                    sigma=error_observ[1:], absolute_sigma=True)
    
    expect_fit = scale_law(length, scale_val[0], scale_val[1])
    fit_plot_val = scale_law(length, copt[0], copt[1])
    
    plt.figure(figsize=[15,8])
    plt.errorbar(length, observ, error_observ, label='simulated observable')
    plt.plot(length, fit_plot_val, label='observable fit')
    plt.plot(length, expect_fit, label='expected scaling law')
    plt.title('Observable with fit values A={} and v={}'.format(np.round(copt[0], decimals=4), 
                                                            np.round(copt[1], decimals=4)), fontsize=20)
    plt.xlabel('Polymer length', fontsize=15)
    plt.ylabel('Observable value', fontsize=15)
    plt.grid()
    plt.legend(loc='best', fontsize='large')
    #plt.savefig("Figures/test1_fit_PERM_end2end.pdf")
    plt.show()
    
    return length, copt, ccov, expect_fit, fit_plot_val

def fit_observ(observ, scale_val=np.array([0.771, 3/4])):
    """Automatically fits the simulated observable and then plot the simulated observable 
    against the fit and the expected behaviour of the observable based on the scaling law.
    
    Parameter
    ---------
    observ : nd.array
        The weighted averages of the observable for different polymer lengths
    error_observ : nd.array
        The error of the weighted averages of th observable for different polymer lengths
    scale_val: nd.array
        The expected values for the scaling factor according to literature. Must be an vector
        with two ellements ([a, b]), with the element relating to the scaling law as:
        <observ(L)> = a * L**(2 * b)
    
    Return
    ------
    length : nd.array
        The length of the polymer from 1 to L
    copt : nd.array
        The fitted values for scale_val with copt[0]=a and copt[1]=b
    ccov : nd.array
        The covariance of the fitted values
    expect_fit : nd.array
        The expected weighted average of the observable based on scale_val
    fit_plot_val: nd.array
        The fit of the weighted average of the observable based on copt
    """
    
    length = np.arange(1, len(observ)+1, 1)
    
    copt, ccov = cv(scale_law, length[1:], observ[1:], p0=scale_val)
    
    expect_fit = scale_law(length, scale_val[0], scale_val[1])
    fit_plot_val = scale_law(length, copt[0], copt[1])
    
    plt.figure(figsize=[15,8])
    plt.plot(length, observ, label='simulated observable')
    plt.plot(length, fit_plot_val, label='observable fit')
    plt.plot(length, expect_fit, label='expected scaling law')
    plt.title('Observable with fit values A={} and v={}'.format(np.round(copt[0], decimals=4), 
                                                            np.round(copt[1], decimals=4)), fontsize=20)
    plt.xlabel('Polymer length', fontsize=15)
    plt.ylabel('Observable value', fontsize=15)
    plt.grid()
    plt.legend(loc='best', fontsize='large')
    #plt.savefig("Figures/test1_fit_PERM_end2end.pdf")
    plt.show()
    
    return length, copt, ccov, expect_fit, fit_plot_val


