using Pkg
using Combinatorics
using LinearAlgebra
using SparseArrays
using Arpack
using Plots
#Pkg.add("Combinatorics")
#Pkg.add("CUDA")
using CUDA

function Base_ind(Nsites)
    a = vec(collect(Base.Iterators.product(Base.Iterators.repeated(1:4, Nsites)...)))
    inds = SparseVector{Int}[]
    for j in a
        push!(inds, [i for i in j])
    end
    return inds
end

function Symmetries(conserve_Sz, conserve_Nt, Sz, Nt, labels)
    labelSym = SparseVector{Int}[]
    if conserve_Nt == true && conserve_Sz == false
        for lab in labels
            Ne_state = 0
            for part in lab
                if (part == 2 || part == 3)
                    Ne_state += 1
                elseif (part == 4) 
                    Ne_state += 2
                end
            end
            if (Ne_state == Nt)
                push!(labelSym, lab)
            end   
        end
    elseif conserve_Nt == false && conserve_Sz == true
        for lab in labels
            Sz_state = 0
            for part in lab
                if (part == 2)
                    Sz_state += 1
                elseif (part == 3) 
                    Sz_state += -1
                end
            end
            if (Sz_state == Sz)
                push!(labelSym, lab)
            end   
        end
    elseif conserve_Nt == true && conserve_Sz == true
        for lab in labels
            Sz_state, Ne_state = 0, 0
            for part in lab
                if (part == 2)
                    Sz_state += 1
                    Ne_state += 1
                elseif (part == 3)
                    Sz_state += -1
                    Ne_state += 1
                elseif (part == 4)
                    Ne_state += 2
                end
            end
            if (Sz_state == Sz && Ne_state == Nt)
                push!(labelSym, lab)
            end   
        end
    end
    return labelSym
end

function C(site, spin, state)
    new_state = copy(state)
    if spin == 2 && state[site] == 2
        new_state[site] = 1
    elseif spin == 3 && state[site] == 3
        new_state[site] = 1
    elseif spin == 2 && state[site] == 4
        new_state[site] = 3
    elseif spin == 3 && state[site] == 4
        new_state[site] = 2
    else
        new_state[site] = 0
    end

    return new_state
end

function Cdag(site, spin, state)
    new_state = copy(state)
    if spin == 2 && state[site] == 1
        new_state[site] = 2
    elseif spin == 3 && state[site] == 1
        new_state[site] = 3
    elseif spin == 2 && state[site] == 3
        new_state[site] = 4
    elseif spin == 3 && state[site] == 2
        new_state[site] = 4
    else
        new_state[site] = 0
    end
    
    return new_state
end

function Hint(Nsites, labels, U)
    Hdim = length(labels)
    UHam = spzeros(Float32, Hdim, Hdim)
    for i = 1:Nsites
        j = 1
        for stateket in labels
            Nup = Cdag(i, 2, C(i, 2, stateket))
            Nupdn = Cdag(i, 3, C(i, 3, Nup))
            k = 1
            for statebra in labels
                if Nupdn == statebra
                    UHam[j, k] += U
                end
                k += 1
            end
            j += 1
        end
    end
    return UHam
end

function Honsite(Nsites, labels, ed)
    Hdim = length(labels)
    edHam = spzeros(Float32, Hdim, Hdim)
    for i = 1:Nsites
        j = 1
        for stateket in labels
            Nup = Cdag(i, 2, C(i, 2, stateket))
            k = 1
            for statebra in labels
                if Nup == statebra
                    edHam[j, k] += -ed
                end
                k += 1
            end
            Ndn = Cdag(i, 3, C(i, 3, stateket))
            k = 1
            for statebra in labels
                if Ndn == statebra
                    edHam[j, k] += -ed
                end
                k += 1
            end
            j += 1
        end
    end
    return edHam
end

function Hhopping(Nsites, labels, t)
    Hdim = length(labels)
    hopHam = spzeros(Float32, Hdim, Hdim)
    for i = 1:(Nsites-1)
        j = 1
        for stateket in labels
            Nup = Cdag(i, 2, C(i+1, 2, stateket))
            k = 1
            for statebra in labels
                if Nup == statebra
                    hopHam[j, k] += -t
                end
                k += 1
            end
            Ndn = Cdag(i, 3, C(i+1, 3, stateket))
            k = 1
            for statebra in labels
                if Ndn == statebra
                    hopHam[j, k] += -t
                end
                k += 1
            end
            j += 1
        end
    end
    for i = 1:(Nsites-1)
        j = 1
        for stateket in labels
            Nup = C(i, 2, Cdag(i+1, 2, stateket))
            k = 1
            for statebra in labels
                if Nup == statebra
                    hopHam[j, k] += -t
                end
                k += 1
            end
            Ndn = C(i, 3, Cdag(i+1, 3, stateket))
            k = 1
            for statebra in labels
                if Ndn == statebra
                    hopHam[j, k] += -t
                end
                k += 1
            end
            j += 1
        end
    end
    return hopHam
end

Nsitios = 4
t = 1
U = 5
mu = U/2

indices = Base_ind(Nsitios)
indices_sim = Symmetries(true, true, Nsitios%2, Nsitios, indices)

Hdim = length(indices_sim)

Hmu = Honsite(Nsitios, indices_sim, mu)
@time Ht = Hhopping(Nsitios, indices_sim, t)
HU = Hint(Nsitios, indices_sim, U)

H = Hmu + Ht + HU + sparse(I, Hdim, Hdim)*U*Nsitios/4