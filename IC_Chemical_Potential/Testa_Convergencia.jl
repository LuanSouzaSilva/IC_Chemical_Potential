using Pkg


Pkg.add("ITensors")
Pkg.add("Plots")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("SparseArrays")
Pkg.add("Arpack")
Pkg.add("LinearAlgebra")
Pkg.add("ProgressBars")



using ITensors, Plots, SparseArrays, Arpack, CSV, DataFrames, LinearAlgebra, ProgressBars


function Hubbard_DMRG(Nsites, t, U, ed, bond_dim)
    sites = siteinds("Electron", Nsites, conserve_nf = false)

    os = OpSum()
    for j=1:Nsites
        os += ed, "Nup",j
        os += ed,"Ndn",j
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

    psi0 = random_mps(sites; linkdims=10) #Esse comando de Random_MPS so funciona se conserve_qns = false

    #state = [isodd(n) ? "Up" : "Dn" for n=1:Nsites]
    #psi0 = MPS(sites,state)

    nsweeps = length(bond_dim)#6
    maxdim = bond_dim#[10, 20, 100, 200, 400, 800]
    cutoff = [1E-12]

    GS_energy, GS = dmrg(H,psi0;nsweeps,maxdim,cutoff, outputlevel = 0)

    Nexp = expect(GS, "Ntot")

    return GS_energy, Nexp
end

function Npart_DMRG(Nsites, t, U, ed, Npart, bond_dim)
    sites = siteinds("Electron", Nsites, conserve_qns = true)

    os = OpSum()
    for j=1:Nsites
        os += ed, "Nup",j
        os += ed,"Ndn",j
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
    psi0 = MPS(sites,state)

    nsweeps = length(bond_dim)#6
    maxdim = bond_dim#[10, 20, 100, 200, 400, 800]
    cutoff = [1E-12]

    GS_energy, GS = dmrg(H,psi0;nsweeps,maxdim,cutoff, outputlevel = 0)

    return GS_energy

end

function Find_EGS(Nsites, Ed_arr, t, U, bond_dim)
    EGSs = zeros(Float64, length(Ed_arr), 3)
    for j in tqdm(eachindex(Ed_arr))
        ed = Ed_arr[j]

        x, y = Hubbard_DMRG(Nsites, t, U, ed, bond_dim)

        local npart_GS = round(sum(y))

        local  EGS_minus = Npart_DMRG(Nsites, t, U, ed, npart_GS - 1, bond_dim)
        local  EGS = Npart_DMRG(Nsites, t, U, ed, npart_GS, bond_dim)
        local EGS_plus = Npart_DMRG(Nsites, t, U, ed, npart_GS + 1, bond_dim)

        #print("\n", x, "  --------------  ", sum(y), "\n")
        #print(EGS_minus, "  --------------  ", EGS_plus, "\n")
        #print(mu, "\n\n")
        EGSs[j, 1] = EGS_minus
        EGSs[j, 2] = EGS
        EGSs[j, 3] = EGS_plus

    end

    return EGSs
end

function Conta_Erro(Nsites, U, Ed_arr, bond_dim)
    @time GSEnergy = Find_EGS(Nsites, Ed_arr, 1, U, bond_dim)

    global contador = 0
    for i in eachindex(Ed_arr)
        if GSEnergy[i, 1] < GSEnergy[i, 2] || GSEnergy[i, 3] < GSEnergy[i, 2]
            global contador += 1
            print(GSEnergy[i, 1], " --------------------- ", GSEnergy[i, 2], "\n")
            print(GSEnergy[i, 3], " --------------------- ", GSEnergy[i, 2], "\n\n")
        end

    end

    print(contador, "\n")
    return contador
end


nsites = 40
u = 10
ed_arr =  LinRange(-0.9*u, -0.7*u, 10) #-u/2 +- u/2

#[10, 20, 100, 200, 400, 800]
bondD = [[10, 20], [10, 20, 100], [10, 20, 100, 200], [10, 20, 100, 200, 400], [10, 20, 100, 200, 400, 800]]

EGS_m1 = zeros(Float64, length(bondD[:, 1]))
EGS_ = zeros(Float64, length(bondD[:, 1]))
EGS_p1 = zeros(Float64, length(bondD[:, 1]))

EGS_apr = zeros(Float64, length(bondD[:, 1]))

for i in tqdm(eachindex(bondD))

    #local count = Conta_Erro(nsites, u, ed_arr, bondD[i])
    EGS_sQN, N_exp = Hubbard_DMRG(nsites, 1, u, -0.8*u, bondD[i])
    npart_GS_apr = round(sum(N_exp))

    EGS = Npart_DMRG(nsites, 1, u, -0.8*u, npart_GS_apr, bondD[i])
    EGSm1 = Npart_DMRG(nsites, 1, u, -0.8*u, npart_GS_apr-1, bondD[i])
    EGSp1 = Npart_DMRG(nsites, 1, u, -0.8*u, npart_GS_apr+1, bondD[i])

    EGS_m1[i] = EGSm1
    EGS_p1[i] = EGSp1
    EGS_[i] = EGS

    EGS_apr[i] = EGS_sQN
end


#p = plot([20, 100, 200, 400, 800], EGS_)
#p = plot!([20, 100, 200, 400, 800], EGS_apr)
#p = plot([20, 100, 200, 400, 800], EGS_m1)
#p = plot!([20, 100, 200, 400, 800], EGS_p1)
#display(p)

df = DataFrame(EGS_without_QN = EGS_apr,
                EGS_with_QN = EGS_, 
                EGS_minus_1 = EGS_m1,
                EGS_plus_1 = EGS_p1)

CSV.write("ConvergenceU10N40.csv", df)

#p = scatter(1/2 .+ ed_arr./u, Chem_Pot./u)
#display(p)