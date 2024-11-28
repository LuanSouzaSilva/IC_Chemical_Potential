using Pkg

Pkg.add("Combinatorics")

using Combinatorics
using LinearAlgebra
using SparseArrays
using Arpack
using Plots
#Pkg.add("CUDA")
#using CUDA

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

function FindNGS(Nsites, labels, ed, t, U)
    E = zeros(Float64, 0)
    for i = 1:(2*Nsites)
        indices_sim = Symmetries(false, true, i%2, i, labels)

        Hdim = length(indices_sim)

        Hed = Honsite(Nsites, indices_sim, ed)
        Ht = Hhopping(Nsites, indices_sim, t)
        HU = Hint(Nsites, indices_sim, U)

        Ham = Ht + HU + Hed + sparse(I, Hdim, Hdim)*U*Nsites/4

        try
            eigvals, _= eigs(Ham, nev = 1, which=:SR, ritzvec = true)
            append!(E, only(eigvals))
        catch
            Ham = Matrix(Ham)
            data = eigmin(Ham)
            #eigvals = data.values
            #eigGS = first(eigvals)
            append!(E, data)
            #print("É, não deu nao!", '\n')
        end
        #print(only(eigvals), '\n')
    end
    _,  ind = findmin(E)

    NGS = ind
        
    return NGS
end

function FindEn(Nsites, labels, ed, t, U, Ne)
    indices_sim = Symmetries(false, true, Ne%2, Ne, labels)

    Hdim = length(indices_sim)

    Hed = Honsite(Nsites, indices_sim, ed)
    Ht = Hhopping(Nsites, indices_sim, t)
    HU = Hint(Nsites, indices_sim, U)

    Ham = Ht + HU + Hed #+ sparse(I, Hdim, Hdim)*U*Nsites/4

    try
        eigvals, _= eigs(Ham, nev = 1, which=:SR, ritzvec = true)
        return only(eigvals)
    catch
        Ham = Matrix(Ham)
        data = eigmin(Ham)
        return data
        #print("É, não deu nao!", '\n')
    end

end

Nsitios = 6

indices = Base_ind(Nsitios)
indices_sim = Symmetries(false, true, Nsitios%2, Nsitios, indices)

Hdim = length(indices_sim)

t = 1
U_arr = [0., 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 10., 12., 14., 16., 18., 20., 30., 40., 50.]  #@time

for i in eachindex(U_arr)
    local E = FindEn(Nsitios, indices, U_arr[i]/2, t, U_arr[i], Nsitios)
    print("\n", E, " --- ", U_arr[i], "\n\n")
end

function rolha()
    onsite = LinRange(0, U, 100)

    #onsite = -U/2

    mu = zeros(Float64, 0)
    os = zeros(Float64, 0)
    contador = 0
    for ed in onsite
        Numgs = FindNGS(Nsitios, indices, ed, t, U)
        Np1 = Numgs + 1
        Nm1 = Numgs - 1

        Ep1 = FindEn(Nsitios, indices, ed, t, U, Np1)
        Em1 = FindEn(Nsitios, indices, ed, t, U, Nm1)

        try
            aux = (Ep1 - Em1)/2
            append!(mu, aux)
            append!(os, ed)
        catch
            print(" ")
        end
        #print(Np1, " - ", Nm1, '\n')
        contador += 1
        print(round(100*contador/length(onsite)), '\n')
    end
    return mu, os
end

#mu, onsite = rolha()

#using CSV

#print(mu)

#x = (1/2 .-onsite./U)
#ind_ = sortperm(x) 
#print(x[ind_])
#CSV.write("ED_N2U5.csv", (mu = mu./U, x))

#print(mu)

#p = plot(x[ind_], .-mu./U, seriestype=:scatter)
#display(p)

#eigvals = data.values
#print(eigvals)

#p = scatter(eigvals)
#display(p)
#p = heatmap(H, yflip = true)
#display(p)