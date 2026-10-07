module Metropolis_MonteCarlo
#include("lattice.jl")
using ..UnitCell_Lattices
using StatsBase
export ising_energy, monte_carlo_simulation, ising_energy_change, kitaev_energy_change, WL_sweep


#Function to calculate the Ising Energy
function ising_energy(lattice::Lattice,J::Vector{Float64},h::Float64) #Energy function of the classical Ising
    E = 0.0
    for site in lattice.sites
        i = site.index
        s_iz = site.spin[3]
        for ngh in lattice.neighbors[i]
            if i < ngh
                s_jz = lattice.sites[ngh].spin[3]
                E -= J[3] * s_iz * s_jz
            end
        end
        E -= h * s_iz
    end
    return E
end

#----------------


# Hamiltonians
function kitaev_energy_change(site::Site,lattice::Lattice,J::Vector{Float64})
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

function ising_energy_change(site::Site,lattice::Lattice,J::Vector{Float64})
    sz = site.spin[3]
    ns = sum(lattice.sites[j].spin[3] for j in lattice.neighbors[site.index])
    de = 2*J[3]*sz*ns
    return de
end
#Single Metropolis flip and Metropolis sweep ~O(N)
function metropolis_step(lattice::Lattice, T::Float64, J::Vector{Float64}, h::Float64, kb=Float64(1.0),model = ising_energy_change)
    i = rand(1:lattice.num_sites)
    site = lattice.sites[i]
    s = site.spin[3]
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
#----------------

# Wang-Landau Sweep

function propose_config(c::Lattice)
    proposed = deepcopy(c)
    i = rand(1:proposed.num_sites)
    proposed.sites[i].spin[3] *= -1
    return proposed
end

function nearest_energy(E_spec::Vector{Float64},E::Float64)
    clossest_indx = findfirst(==(E), E_spec)
    isnothing(clossest_indx) && error("Computed energy $E is missing from E_spec")
    clossest_E = E_spec[clossest_indx]
    return  clossest_indx, clossest_E
end

function isFlat(h::Vector{Int},flatness::Float64 = 0.8)
    h_bar = mean(h)
    h_min = minimum(h)
    if h_min > flatness*h_bar
        return true
    else
        return false
    end
end

function WL_sweep(c0::Lattice,E_spec::Vector{Float64},J::Vector{Float64},h::Float64)
    log_g = zeros(length(E_spec))
    H = zeros(Int,length(E_spec))
    log_f = 1.0
    current_config = c0
    while log_f > 1e-3
        n_config = propose_config(current_config)
        E_config_current = ising_energy(current_config,J,h)
        E_config_proposed = ising_energy(n_config,J,h)
        Eindx_current, E_current = nearest_energy(E_spec,E_config_current)
        Eindx_proposed, E_proposed = nearest_energy(E_spec,E_config_proposed)
        g_current = log_g[Eindx_current]
        g_proposed = log_g[Eindx_proposed]
        prob_accept = min(0.0,g_current-g_proposed)
        if log(rand()) < prob_accept
            current_config = n_config
            occupied_indx = Eindx_proposed
        else
            occupied_indx = Eindx_current
        end
        
        log_g[occupied_indx] += log_f
        H[occupied_indx] += 1

        #Check if the histogram is flat
        if isFlat(H)
            #DO measurements here
            println("Flat histogram is reached!!")
            H = zeros(Int,length(E_spec))
            log_f /=2
            println(log_f)
        end
        
    end

    #Normalize log_g
    log_g_norm = log_g .+ log(2) .- log_g[1] #THe first element corresponds to the log_g(E_min)
    return log_g
end



#-----------

#Simulation Functions

function WL_simulate(lattice::Lattice, J::Vector{Float64}, h::Float64, T::Float64, num_steps_thermalization::Int)

end


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