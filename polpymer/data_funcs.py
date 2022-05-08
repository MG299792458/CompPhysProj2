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

