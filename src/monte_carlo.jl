module Metropolis_MonteCarlo

using ..UnitCell_Lattices
using StatsBase
export ising_energy, monte_carlo_simulation


#Function to calculate the Ising Energy
function ising_energy(lattice::Lattice,J::Vector{Float64},h::Float64) #Energy function of the classical Ising
    E = 0.0
    for site in lattice.sites
        i = site.index
        s_iz = site.spin[3]
        for ngh = lattice.neighbors[i]
            s_jz = lattice.sites[ngh].spin[3]
            E -= J[1] * s_iz * s_jz
        end
        E -= h * s_iz
    end
    return E
end

# Hamiltonians
function kitaev_energy(site::Site,lattice::Lattice,J::Vector{Float64})
    sx = site.spin[1]
    ns_x = sum(lattice.sites[j].spin[1] for j in lattice.neighbors[site.index])
    de_x = 2*J[1]*sx*ns_x

    sy = site.spin[2]
    ns_y = sum(lattice.sites[j].spin[2] for j in lattice.neighbors[site.index])
    de_y = 2*J[2]*sy*ns_y

    sz = site.spin[3]
    ns_z = sum(lattice.sites[j].spin[3] for j in lattice.neighbors[site.index])
    de_z = 2*J[3]*sz*ns_z

    return de_x + de_y + de_z
end
#Single Metropolis flip and Metropolis sweep ~O(N)
function metropolis_step(lattice::Lattice, T::Float64, J::Vector{Float64}, h::Float64, kb=Float64(1.0),model = kitaev_model)
    rand_i = rand(1:lattice.num_sites)
    i = rand(1:lattice.num_sites)
    site = lattice.sites[i]
    dE = model(site,lattice,J)
    if dE <= 0 || rand() < exp(-dE*kb / T)
        lattice.sites[i].spin[3] = -s
    end
end

function metropolis_sweep!(lattice::Lattice, T::Float64, J::Vector{Float64}, h::Float64, kb=Float64(1.0))
    for _ in 1:lattice.num_sites
        metropolis_step(lattice, T,J,h, kb)
    end
end

#Simulation Function
function monte_carlo_simulation(lattice::Lattice, T::Float64, J::Vector{Float64}, h::Float64, num_steps::Int, num_equilibration::Int,num_data_skip::Int,kb=Float64(1.0))
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