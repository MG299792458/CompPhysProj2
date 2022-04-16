""" Polpymer Core Functions

This script contains all the functions used to simulate the creation of random
polymers by using a Monte-Carlo method.

Results of the simulation can be analysed using the functions in data_funcs.py
"""


# Module imports
from tkinter import E
import matplotlib.pyplot as plt
from polpymer.core_funcs import Polymer, Monomer
import random
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


def grow_polymer(dims, origin, L: int):
        """randomly grows a polymer up to a length of L or until it can't grow anymore and stores the number of growth option for each growth step to determine the weigth of the polymer
        
        """
        first_monomer = Monomer(random.choice(range(4)))
        polymer = Polymer((10,10), (5,5), first_monomer)
        m = 4*np.ones(L)
        for i in range(1,L):
            if m[i-1] > 0:
                grow_directions = [0,1,2,3]
                grow_options = [0,1,2,3]
                for j in grow_directions:
                    proposed_monomer = Monomer(j)
                    proposed_monomer.location = polymer.chain_end
                    proposed_monomer.calculate_end()
                    if polymer.conflict(proposed_monomer):
                        m[i] -= 1
                        grow_options.remove(j)
                if len(grow_options) > 0:
                    polymer.add_monomer(random.choice(grow_options))
            else:
                print("The polymer couldn't grow to length L")
                break
        return m, polymer          
                
                
            
        
        