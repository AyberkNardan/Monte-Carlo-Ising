# A Module to do Montecarlo Simulation of Classical Ising Model
Implementation of Metropolis step Monte Carlo simulations for classical spin systems. The module works on arbitrary lattices. Further models, like the Heisenberg Model, can be implemented easily by edding the energy calculation into the `metropolis_ising.jl`. 

One can also use the code just for constructing any type of lattice with different spin models. For now, the lattice construction only works for the PBC.
