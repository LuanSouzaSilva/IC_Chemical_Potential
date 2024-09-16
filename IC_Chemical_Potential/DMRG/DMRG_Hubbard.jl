using Pkg

Pkg.add("ITensors")
Pkg.add("Plots")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("SparseArrays")
Pkg.add("Arpack")
Pkg.add("LinearAlgebra")
Pkg.add("ProgressBars")

Pkg.instantiate()

using ITensors, Plots, SparseArrays, Arpack, CSV, DataFrames, LinearAlgebra, ProgressBars


function Hubbard_DMRG(Nsites, t, U, ed)
    sites = siteinds("Electron", Nsites, conserve_nf = false)

    os = OpSum()
    for j=1:Nsites
        os += -ed, "Nup",j
        os += -ed,"Ndn",j
    end

    for j=1:Nsites-1
        os += -t, "Cdagup",j+1, "Cup", j
        os += -t, "Cdagup",j, "Cup", j+1

        os += -t,"Cdagdn",j+1, "Cdn", j
        os += -t,"Cdagdn",j, "Cdn", j+1
    end

    for j=1:Nsites
        os += U, "Nupdn",j
    end

    H = MPO(os,sites)

    psi0 = randomMPS(sites;linkdims=10) #Esse comando de Random_MPS so funciona se conserve_qns = false

    #state = [isodd(n) ? "Up" : "Dn" for n=1:Nsites]
    #psi0 = MPS(sites,state)

    nsweeps = 6
    maxdim = [10, 20, 100, 200, 400, 800]
    cutoff = [1E-12]

    GS_energy, GS = dmrg(H,psi0;nsweeps,maxdim,cutoff, outputlevel = 0)

    Nexp = expect(GS, "Ntot")

    return GS_energy, Nexp
end

function Npart_DMRG(Nsites, t, U, ed, Npart)
    sites = siteinds("Electron", Nsites, conserve_qns = true)

    os = OpSum()
    for j=1:Nsites
        os += -ed, "Nup",j
        os += -ed,"Ndn",j
    end

    for j=1:Nsites-1
        os += -t, "Cdagup",j+1, "Cup", j
        os += -t, "Cdagup",j, "Cup", j+1

        os += -t,"Cdagdn",j+1, "Cdn", j
        os += -t,"Cdagdn",j, "Cdn", j+1
    end

    for j=1:Nsites
        os += U, "Nupdn",j
    end

    H = MPO(os,sites)

    state = ["Emp" for n in 1:Nsites]
    p = Npart
    for i in Nsites:-1:1
        if p > i
            state[i] = "UpDn"
            p -= 2
        elseif p > 0
            state[i] = (isodd(i) ? "Up" : "Dn")
            p -= 1
        end
    end
    psi0 = randomMPS(sites,state, linkdims = 10)

    nsweeps = 6
    maxdim = [10, 20, 100, 200, 400, 800]
    cutoff = [1E-12]

    GS_energy, GS = dmrg(H,psi0;nsweeps,maxdim,cutoff, outputlevel = 0)

    return GS_energy

end

function Chemical_Potential(Nsites, Ed_arr, t, U)
    Chem_Pot = zeros(Float64, 0)
    for ed in tqdm(Ed_arr)

        x, y = Hubbard_DMRG(Nsites, t, U, ed)

        local npart_GS = round(sum(y))

        local  EGS_minus = Npart_DMRG(Nsites, t, U, ed, npart_GS - 1)
        local EGS_plus = Npart_DMRG(Nsites, t, U, ed, npart_GS + 1)

        local mu = (EGS_plus-EGS_minus)/2

        #print("\n", x, "  --------------  ", sum(y), "\n")
        #print(EGS_minus, "  --------------  ", EGS_plus, "\n")
        #print(mu, "\n\n")
        append!(Chem_Pot, mu)

    end

    return Chem_Pot
end

function Roda_N(U, nsites_arr, filenames)
    ed_arr =  LinRange(0.9*U, 0.75*U, 100) #u/2 +- u/2

    for i in eachindex(nsites_arr)
        Chem_Pot = Chemical_Potential(nsites_arr[i], ed_arr, 1, U)

        df = DataFrame(Onsite_Energy = ed_arr,
                    Chemical_Potential = Chem_Pot)
        
        CSV.write(filenames[i], df)

    end
end

Nsites = [10, 20, 40, 60, 100]
Filenames = ["U10N10.csv", "U10N20.csv", "U10N40.csv", "U10N60.csv", "U10N100.csv"]

u = 10

Roda_N(u, Nsites, Filenames)