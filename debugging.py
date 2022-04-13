from polpymer.data_funcs import Polymer, Monomer

start_monomer = Monomer(0)

polymer = Polymer((10,10), (0,0), start_monomer)
print(polymer)
polymer.add_monomer(3)
print(polymer)