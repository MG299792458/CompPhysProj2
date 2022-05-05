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


def grow_polymer(dims, origin, L: int):
        """randomly grows a polymer up to a length of L or until it can't grow anymore and stores the number of growth option for each growth step to determine the weigth of the polymer

        """
        polymer = Polymer(dims, origin)
        m = np.zeros(L-1)
        grow_directions = [0,1,2,3]
        for i in range(L-1):
            grow_options = [0,1,2,3]
            for j in grow_directions:
                proposed_monomer = Monomer(j)
                proposed_monomer.location = polymer.chain_end
                proposed_monomer.calculate_end()
                if polymer.conflict(proposed_monomer):
                    grow_options.remove(j)
            if len(grow_options) > 0:
                    m[i] = len(grow_options)
                    polymer.add_monomer(choice(grow_options))
            else:
                print("The polymer grew to length {}".format(i+1))
                break
        m = m[m!=0]
        return m, polymer


def find_polymer(dims,
                 origin,
                 L):
    """find a polymer that has the desired lenght L

    Parameters
    ----------
    dims : Tuple[int, int]
        Amount of nodes in either x and y direction, unused for now
    origin : Tuple[int,int]
        Starting node of the first monomer
    L : Tuple[int]
        Length of each polymer

    Return
    ------
    m : nd.array
        the weight of the polymer
    polymer: polymer of length L
    """

    n = 0

    while n != L:
        m, polymer = grow_polymer(dims, origin, L)
        n = polymer.chain_length

    return m, polymer


def read_polymer(polymer: object, m):
    """Reading location of the nodes of a single polymer and the corresponding weight.

    Parameters
    ----------
    polymer :
    m : nd.array
        the number of unoccupied lattice sites at each polymer node

    Return
    ------
    x_ : nd.array
        the x coordinate for each polymer node
    y_ : nd.array
        the y coordinate for each polymer node
    w_ : nd.array
        the weight of each polymer node
    """

    L = polymer.chain_length #number of nodes
    x_ = np.array([])
    y_ = np.array([])
    w_ = np.array([])

    difference = (polymer[0].location[0]-polymer[-1].location[0], polymer[0].location[1]-polymer[-1].location[1])
    end_to_end = (difference[0]**2 + difference[1]**2)

    for monomer in polymer:
        start = monomer.location
        end = monomer.end_location

        x_ = np.append(x_, start[0])
        y_ = np.append(y_, start[1])

    for i in range(L):
        w_ = np.append(w_, np.prod(m[0:i+1]))

    return x_, y_, w_


def observ_polymer(x_, y_):
    """Calculating the end_to_end distance and the the radius of gyration of a single polymer

    Parameter
    ---------
    x_ : nd.array
        the x coordinate for each polymer node
    y_ : nd.array
        the y coordinate for each polymer node

    Return
    ------
    end_to_end : nd.array
        the end_to_end distance for all the polymer nodes
    gyration : nd.array
        the radius of gyration for all the polymer nodes
    """

    L = np.size(x_) #number of nodes

    end_to_end = np.array([])
    gyration = np.array([])

    for i in range(L-1):
        end_to_end_i = (x_[0] - x_[i+1])**2 + (y_[0] - y_[i+1])**2
        end_to_end = np.append(end_to_end, end_to_end_i)

        cm_x = 1/(i+1) * np.sum(x_[0:i+1])
        cm_y = 1/(i+1) * np.sum(y_[0:i+1])

        gyration_x_i = 1/(i+1) * np.sum((x_[0:i+1] - cm_x)**2)
        gyration_y_i = 1/(i+1) * np.sum((y_[0:i+1] - cm_y)**2)
        gyration_i = gyration_x_i + gyration_y_i
        gyration = np.append(gyration, gyration_i)

    return end_to_end, gyration


def generate_N_polymers(N: int, L: int, dim, origin):
    """fuction to generate N polymers of length L

    Parameter
    ---------
    N : int
        number of polymers to generate
    L : int
        length of the polymers generated

    Return
    ------
    end_to_end : nd.array
        N x L matrix where the (i,j) element represents the end_to_end distance of polymer i with length j+1
    gyration : nd.array
        N x L matrix where the (i,j) element represents the radius of gyration of polymer i with length j+1
    w : nd.array
        N x L matrix where the (i,j) element represents the weight of polymer
    """

    end_to_end = np.zeros((N,L-1))
    gyration = np.zeros((N,L-1))
    w = np.zeros((N,L-1))

    for i in range(N):
        m, polymer = find_polymer(dim, origin, L)

        x_i, y_i, w_i = read_polymer(polymer, m)
        end_to_end_i, gyration_i = observ_polymer(x_i, y_i)

        end_to_end[i,:] = end_to_end_i
        gyration[i,:] = gyration_i
        w[i,:] = w_i[0:-1]

    return end_to_end, gyration, w


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

