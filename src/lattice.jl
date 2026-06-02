module UnitCell_Lattices

export UnitCell, BondRule, Site, Lattice, build_lattice

struct UnitCell
    lattice_vecs::Matrix{Float64} # Matrix where each column is a lattice vector
    basis_vecs::Matrix{Float64} # Matrix where each column is a basis vector
end

struct BondRule
    site1::Int # Index of the first site in the current unit cell
    site2::Int # Index of the second site in the target unit cell
    bond_vec::Vector{Int} # Vector connects two neighboring unit cells 
end

mutable struct Site
    index::Int
    cell::Vector{Int}
    basis::Int
    position::Vector{Float64}
    spin::Vector{Float64}
end

mutable struct Lattice
    unit_cell::UnitCell
    bond_rules::Vector{BondRule}
    lattice_size::Vector{Int} # Number of unit cells in each direction [Lx, Ly, Lz]
    sites::Vector{Site}
    neighbors::Vector{Vector{Int}}
    num_sites::Int
    num_basis::Int
end

# I want to reach each site in the lattice with a unique index. Thus, I will define a function
# that takes the cell [nx,ny,nz], corresponding basis index and the lattice size [Lx, Ly, Lz] as well as
# To do that, I will also define an indexing function for the cell [nx,ny,nz] as
# cell_index = 1 + nx + Lx*ny + Lx*Ly*nz + ... (if we have higher dimensions)
# Then, the site index becomes site_index = (cell_index-1) * num_basis + basisi_index
# Here, we substract 1 from cell_index to make the first site index start from 1 as julia requires
# (1-1)*num_basis + 1 = 1

function cell_index(cell::Vector{Int}, lattice_size::Vector{Int})
    index = 1
    step = 1
    for i in 1:length(cell)
        index += cell[i]*step
        step *= lattice_size[i]
    end
    return index
end

function site_index(cell::Vector{Int}, basis::Int, lattice_size::Vector{Int}, num_basis::Int)
    cellidx = cell_index(cell, lattice_size)
    return (cellidx - 1) * num_basis + basis
end


# Now we can build the lattice

function build_lattice(
    unitcell::UnitCell,
    bond_rules::Vector{BondRule},
    lattice_size::Vector{Int},
    initial_spin = [0.0,0.0,0.0], #Defining the default spin as a zero 3D vector.
    boundary_condition = "periodic" 
)

    dim = length(lattice_size)
    num_basis = size(unitcell.basis_vecs,2)
    num_cell = prod(lattice_size)
    num_sites = num_cell * num_basis

    sites = Vector{Site}(undef, num_sites) #Creating the list of sites in the lattice with number of sites
    neighbors = [Int[] for i in 1:num_sites] #Creating the neighbors list for each site

    ranges = [0:(lattice_size[i]-1) for i in 1:dim] #Creating the range of cell indices for each dimension
    # This ranges will return cartesian product of cell indices via the Iterators.product function

    # creating the sites in the lattice
    for cell in Iterators.product(ranges...)
        cell = [cell[i] for i in 1:length(cell)] #Converting the cell tuple to a vector because Iterators.product returns a tuple

        for b in 1:num_basis
            b_vec = unitcell.basis_vecs[:,b]
            Rn = unitcell.lattice_vecs * cell
            site_pos = Rn + b_vec
            siteidx = site_index(cell,b,lattice_size,num_basis)
            sites[siteidx] = Site(siteidx, cell, b, site_pos, initial_spin)
        end
    end

    # defining the neighbors according to the bond rules
    for cell in Iterators.product(ranges...)
        cell = [cell[i] for i in 1:length(cell)]

        for rule in bond_rules
            initial_cell = cell
            target_cell = cell .+ rule.bond_vec
            
            #by the periodic boundary condition,
            if boundary_condition == "periodic" #for now it only supports periodic boundary condition
                target_cell = mod.(target_cell, lattice_size) # returns mod(target_cell[i], lattice_size[i]) for each i
            end
            i = site_index(initial_cell, rule.site1, lattice_size, num_basis)
            j = site_index(target_cell, rule.site2, lattice_size, num_basis)
            push!(neighbors[i], j)
            push!(neighbors[j], i) # By the lattice symmetry.
        end
    end

    return Lattice(unitcell, bond_rules, lattice_size, sites, neighbors, num_sites, num_basis)

end


end

