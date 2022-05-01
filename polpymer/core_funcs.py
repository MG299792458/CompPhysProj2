"""Polpymer Data Processing

This script contains all the functions used to analyse the physicial quantities
and observables of the simulated polymer.

These functions are written to be used on the results obtained by use of the
functions in core_funcs.py
"""


# Module imports
from typing import Tuple
from random import choice
import numpy as np

# Defining global variables
ANGLE_TO_ADD: list[Tuple[int,int]] = [
    (1,0),
    (0,1),
    (-1,0),
    (0,-1)
]


# Top-level functions and classes
class Monomer:
    """ Single Monomere element part of longe Polymer chain.
    single monomer groups together starting and ending point.
    """
    def __init__(self, ang: int):
        """Initialises Monomer class, takes single argument ang(le).

        Parameters
        ----------
        ang : int
            integer angle argument where 1 correpsonds to 90 degrees, options are
            0, 1, 2, 3 for x, y, -x, -y respectively
        """
        self.angle = ang
        self.location: Tuple[int,int] = None
        self.end_location: Tuple[int,int] = None

    def __str__(self):
        """ Prints description of monomer by specifying start and end points.

        Returns
        -------
        _type_
            _description_
        """
        string: str = "Monomer from {} to {}".format(self.location,
        self.end_location)
        return string

    def calculate_end(self):
        """ Calculates end coordinates of monomer link based on angle and
        starting position

        Raises
        ------
        ValueError
            If function is called without specifying the starting location first
        """
        if self.location is None:
            raise ValueError("Location of end not possible when location is None")
        else:
            add = ANGLE_TO_ADD[self.angle]
            start_loc = self.location
            self.end_location = (start_loc[0]+add[0], start_loc[1]+add[1])

    def calculate_cm(self) -> Tuple[float,float]:
        """Calculates the centre of mass of this monomer in global coordinates

        Returns
        -------
        Tuple[float,float]
            x, y coordinate pair of the centre of mass
        """
        if self.end_location is None:
            self.calculate_end()

        cm: Tuple[float,float] = None
        xcm: float = self.location[0] + (self.end_location[0]-self.location[0])/2
        ycm: float = self.location[1] + (self.end_location[1]-self.location[1])/2

        cm = (xcm, ycm)
        self.mass_centre = cm

        return cm

class Polymer:
    """ Polymer object encapsulates dictionary of monomer objects

    Returns
    -------
    Polymer
        Object containing grouped information on the polymer chain

    Yields
    ------
    Iterable
        Allows for iterating over the monomers in the polymer

    Raises
    ------
    ValueError
        If invalid parameter is passed to the add_monomer function
    Exception
        If adding the monomer creates a self crossing
    """

    chain_length: int = 0
    chain_start: Tuple[int,int] = None
    chain_end: Tuple[int,int] = None

    monomers: dict[int, Monomer] = {}
    claimed_sites: list[Tuple[int,int]] = []

    nodes_locsx: np.ndarray = None
    nodes_locsy: np.ndarray = None
    node_m_vals: list[int,...] = []
    node_weights: list[int,...] = []

    end_to_end: np.ndarray = None
    gyration: np.ndarray = None
    
    pruned: bool = None

    def __init__(self,
        dims: Tuple[int, int],
        origin: Tuple[int,int]):
        """ Initialises the Polymer class with an initial Monomer

        Parameters
        ----------
        dims : Tuple[int, int]
            Amount of nodes in either x and y direction, unused for now
        origin : Tuple[int,int]
            Starting node of the first monomer
        starting_monomer : Monomer
            Monomer object that will start the chain
        """

        self.dimensions = dims
        self.origin = origin

        starting_monomer = Monomer(choice([0,1,2,3]))

        starting_monomer.location = origin
        starting_monomer.calculate_end()

        self.chain_start = origin
        self.chain_end = starting_monomer.end_location

        self.monomers = {}
        self.monomers['monomer_0'] = starting_monomer
        self.chain_length = 1

        self.claimed_sites = [origin, self.chain_end]
        
        self.pruned = False

    def __iter__(self) -> Monomer:
        for i in range(len(self)):
            string: str = 'monomer_{}'.format(str(i))
            yield self.monomers[string]

    def __str__(self):
        string: str = "Polymer chain consisting of {} monomers".format(self.chain_length)
        return string

    def __getitem__(self, item):
        item_str: str = "monomer_{}".format(item+1)
        return self.monomers[item_str]

    def __len__(self):
        return self.chain_length

    def add_monomer(self, ang: int, loc: str = 'end'):
        """ Adds a monomer to the polymer chain

        Parameters
        ----------
        ang : int
            Angle of new monomer with respect to the global orientation
        loc : str, optional
            wether to add the monomer to the starting node or ending node.
            specify with either 'start' or 'end', by default 'end'

        Raises
        ------
        ValueError
            If string loc not 'start' or 'end'
        Exception
            If addition of monomer would create a self-crossing.
        """
        if loc == 'start':
            start_loc = self.chain_start
        elif loc == 'end':
            start_loc = self.chain_end
        else:
            raise ValueError("string location either 'start' or 'end',\
                 default is 'end'")

        proposed_monomer = Monomer(ang)
        proposed_monomer.location = start_loc
        proposed_monomer.calculate_end()

        if not self.conflict(proposed_monomer):
            self.chain_length += 1
            self.monomers['monomer_'+str(self.chain_length-1)] = proposed_monomer
            if start_loc == self.chain_start:
                self.chain_start = proposed_monomer.end_location
            elif start_loc == self.chain_end:
                self.chain_end = proposed_monomer.end_location
            self.claimed_sites.append(self.chain_end)
        else:
            raise Exception("Proposed monomer's end location already a node of polymer")


    def conflict(self, prop_monomer: Monomer) -> bool:
        start = prop_monomer.location
        end = prop_monomer.end_location
        cross: bool = bool(end in self.claimed_sites)
        attach_end: bool = not bool((start == self.chain_start) or (start == self.chain_end))
        close: bool = bool(end == self.chain_start)

        is_conflicting: bool = bool(cross or attach_end or close)
        return is_conflicting

    def observables(self):
        if (self.nodes_locsx is None) or (self.nodes_locsy is None):
            self.compute_node_weights()

        L = self.chain_length
        x_ = self.nodes_locsx
        y_ = self.nodes_locsy

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


    def grow_polymer(self, length):
        """randomly grows a polymer up to a length of L or until it can't grow anymore and stores the number of growth option for each growth step to determine the weigth of the polymer

        """

        m = [4]
        grow_directions = [0,1,2,3]
        for i in range(length-1):
            grow_options = [0,1,2,3]
            for j in grow_directions:
                proposed_monomer = Monomer(j)
                proposed_monomer.location = self.chain_end
                proposed_monomer.calculate_end()
                if self.conflict(proposed_monomer):
                    grow_options.remove(j)
            if len(grow_options) > 0:
                    #m[i] = len(grow_options)
                    m.append(len(grow_options))
                    self.add_monomer(choice(grow_options))
            else:
                print("The polymer grew to length {}".format(i+1))
                break
        #m = m[m!=0]
        self.node_m_vals = m

    def compute_node_weights(self):
        m = self.node_m_vals
        length = self.chain_length

        x_ = np.array([])
        y_ = np.array([])
        w_ = np.array([])

        polymer = self

        difference = (polymer[0].location[0]-polymer[-1].location[0], polymer[0].location[1]-polymer[-1].location[1])
        end_to_end = (difference[0]**2 + difference[1]**2)

        for monomer in polymer:
            start = monomer.location

            x_ = np.append(x_, start[0])
            y_ = np.append(y_, start[1])

        for i in range(length):
            w_ = np.append(w_, np.prod(m[0:i+1]))

        self.nodes_locsx = x_
        self.nodes_locsy = y_
        self.node_weights = w_


class Dish: #As in a Petri-dish

    polymers: list[object,...] = []
    end_to_end: np.ndarray = None
    gyration: np.ndarray = None
    weights: np.ndarray = None

    def __init__(self, dims: Tuple[int,int], origin: Tuple[int,int]):
        self.dimension = dims
        self.origin = origin

    def find_polymer(self, length: int):
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

        dims = self.dimension
        origin = self.origin

        n = 0
        while n != length:
            trial_polymer = Polymer(dims, origin)
            trial_polymer.grow_polymer(length)
            n = trial_polymer.chain_length

        m = trial_polymer.node_m_vals
        return m, trial_polymer

    def generate_N_polymers(self, N: int, length: int):
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

        end_to_end = np.zeros((N,length-1))
        gyration = np.zeros((N,length-1))
        w = np.zeros((N,length-1))

        for i in range(N):
            m, polymer = self.find_polymer(length)

            self.polymers.append(polymer)

            polymer.compute_node_weights()
            x_i = polymer.nodes_locsx
            y_i = polymer.nodes_locsy
            w_i = polymer.node_weights
            end_to_end_i, gyration_i = polymer.observables()

            end_to_end[i,:] = end_to_end_i
            gyration[i,:] = gyration_i
            w[i,:] = w_i[0:-1]

        self.end_to_end = end_to_end
        self.gyration = gyration
        self.weights = w

        return end_to_end, gyration, w
    
    def PERM(self, N: int, cplus: float, L: int):
        cminus = cplus/10
        
        for i in range(N):
            m, polymer = self.find_polymer(1)
            self.polymers.append(polymer)
            
        grow_directions = [0,1,2,3]
            
        for i in range(L-1):
            w = []
            N_polymers = 0
            for polymer in self.polymers:
                m = polymer.node_m_vals
                grow_options = [0,1,2,3]
                if not polymer.pruned:
                    for j in grow_directions:
                        proposed_monomer = Monomer(j)
                        proposed_monomer.location = polymer.chain_end
                        proposed_monomer.calculate_end()
                        if polymer.conflict(proposed_monomer):
                            grow_options.remove(j)
                    if len(grow_options) > 0:
                        m.append(len(grow_options))
                        polymer.add_monomer(choice(grow_options))
                    else:
                        m.append(0)
                    N_polymers += 1
                else:
                    m.append(0)
                polymer.node_m_vals = m
                polymer.compute_node_weights()
                w.append(polymer.node_weights[-1])
                
                
            W_tilde = sum(w)/N_polymers
            W_plus = cplus*W_tilde
            W_minus = cminus*W_tilde
            
            copied_polymers = []
            
            for polymer in self.polymers:
                if polymer.node_weights[-1] < W_minus:
                    if choice([0,1]) == 0:
                        polymer.node_m_vals[-1] = 0
                        polymer.pruned = True
                    else:
                        polymer.node_m_vals[-1] = 2*polymer.node_m_vals[-1]
                    polymer.compute_node_weights()
                elif polymer.node_weights[-1] > W_plus:
                    polymer.node_m_vals[-1] = 0.5*polymer.node_m_vals[-1]
                    polymer.compute_node_weights()
                    copied_polymers.append(polymer)
                    
            for polymer in copied_polymers:
                self.polymers.append(polymer)
                
                
# Checking if the file is ran by itself or imported:
if __name__ == "__main__":
    print("import module with 'from polpymer.data_funcs import *' \
    , instead of running directly")


