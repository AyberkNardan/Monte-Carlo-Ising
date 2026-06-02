module Metropolis_MonteCarlo

using ..UnitCell_Lattices
using StatsBase
export ising_energy, monte_carlo_simulation


#Function to calculate the Ising Energy
function ising_energy(lattice::Lattice,J::Float64,h::Float64) #Energy function of the classical Ising
    E = 0.0
    for site in lattice.sites
        i = site.index
        s_iz = site.spin[3]
        for ngh = lattice.neighbors[i]
            s_jz = lattice.sites[ngh].spin[3]
            E -= J * s_iz * s_jz
        end
        E -= h * s_iz
    end
    return E
end

#Single Metropolis flip and Metropolis sweep ~O(N)
function metropolis_step(lattice::Lattice, T::Float64, J::Float64, h::Float64, kb=Float64(1.0))
    rand_i = rand(1:lattice.num_sites)
    i = rand(1:lattice.num_sites)
    s = lattice.sites[i].spin[3]
    neighbor_sum = sum(lattice.sites[j].spin[3] for j in lattice.neighbors[i])
    dE = 2 * J * s * neighbor_sum + 2 * h * s
    if dE <= 0 || rand() < exp(-dE*kb / T)
        lattice.sites[i].spin[3] = -s
    end
end

function metropolis_sweep!(lattice::Lattice, T::Float64, J::Float64, h::Float64, kb=Float64(1.0))
    for _ in 1:lattice.num_sites
        metropolis_step(lattice, T,J,h, kb)
    end
end

#Simulation Function
function monte_carlo_simulation(lattice::Lattice, T::Float64, J::Float64, h::Float64, num_steps::Int, num_equilibration::Int,num_data_skip::Int,kb=Float64(1.0))
    energies = Float64[]
    magnetizations = Float64[]
    magnetization_squared = Float64[]
    sweep = 0
    for step in 1:num_steps
        for _ in 1:lattice.num_sites
            metropolis_step(lattice, T, J, h,kb)
        end
    sweep += 1
        if step > num_equilibration && step % num_data_skip == 0
            push!(energies, ising_energy(lattice,J,h))
            push!(magnetizations, mean([site.spin[3] for site in lattice.sites]))
            push!(magnetization_squared,mean([site.spin[3] for site in lattice.sites])^2)
        end
    end
    return energies, magnetizations, magnetization_squared
end


end