using Pkg

Pkg.add("ITensors")
Pkg.add("Plots")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("SparseArrays")
Pkg.add("Arpack")
Pkg.add("LinearAlgebra")
Pkg.add("ProgressBars")

include("DMRG_Hubbard.jl")

Pkg.instantiate()

using ITensors, Plots, SparseArrays, Arpack, CSV, DataFrames, LinearAlgebra, ProgressBars

Threads.nthreads() = 64

function Track_N(nsites, ed_arr, t, U)
    Chem_pot = zeros(Float64, 0)
    E1m2_arr = zeros(Float64, 0)
    E1p2_arr = zeros(Float64, 0)
    Npart_arr = zeros(Float64, 0)
    E1_arr = zeros(Float64, 0)


    GS_trial, Npart_trial_exp = Hubbard_DMRG(nsites, t, U, ed_arr[1])
    Npart_trial = round(sum(Npart_trial_exp))

    global GS_Nconserve = Npart_DMRG(nsites, t, U, ed_arr[1], Npart_trial)
    global GS_Nconservep1 = Npart_DMRG(nsites, t, U, ed_arr[1], Npart_trial + 1)
    global GS_Nconservem1 = Npart_DMRG(nsites, t, U, ed_arr[1], Npart_trial - 1)
        
    if GS_Nconserve < GS_Nconservem1 || GS_Nconserve < GS_Nconservep1
        global Npart = Npart_trial
    elseif GS_Nconservem1 < GS_Nconserve
        global Npart = Npart_trial - 1
    else
        global Npart = Npart_trial + 1
    end

    local E1 = Npart_DMRG(nsites, t, U, ed_arr[1], Npart_trial)

    GSp1 = Npart_DMRG(nsites, t, U, ed_arr[1], Npart+1)
    GSm1 = Npart_DMRG(nsites, t, U, ed_arr[1], Npart-1)

    GSp2 = Npart_DMRG(nsites, t, U, ed_arr[1], Npart+2)
    GSm2 = Npart_DMRG(nsites, t, U, ed_arr[1], Npart-2)

    mu = (GSp1 - GSm1)/2

    append!(E1m2_arr, GSm2)
    append!(E1p2_arr, GSp2)

    append!(Chem_pot, mu)

    append!(E1_arr, GS_Nconserve)
    append!(Npart_arr, Npart)
    for  i in tqdm(eachindex(ed_arr)[2:length(ed_arr)])
        local GS_Nconserve = Npart_DMRG(nsites, t, U, ed_arr[i], Npart)
        local GS_Nconservep1 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart + 1)
        local GS_Nconservem1 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart - 1)

        if GS_Nconserve < GS_Nconservem1 && GS_Nconserve < GS_Nconservep1
            Npart += 0
            local E1 = GS_Nconserve
            local GSp1 = GS_Nconservep1
            local GSm1 = GS_Nconservem1
        elseif GS_Nconservem1 < GS_Nconserve
            Npart -= 1
            local E1 = GS_Nconservem1
            local GSp1 = GS_Nconserve
            local GSm1 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart-1)
        elseif GS_Nconservep1 < GS_Nconserve
            Npart += 1
            local E1 = GS_Nconservep1
            local GSm1 = GS_Nconserve
            local GSp1 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart+1)
        end

        GSp2 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart+2)
        GSm2 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart-2)

        mu = (GSp1 - GSm1)/2

        append!(E1m2_arr, GSm2)
        append!(E1p2_arr, GSp2)
        #mu2 = (GSp12 - GSm12)/2

        append!(Chem_pot, mu)

        append!(E1_arr, E1)
        append!(Npart_arr, Npart)

    end

        return Chem_pot, Npart_arr, E1_arr, E1m2_arr, E1p2_arr

end

function Roda_N(U, nsites_arr, filenames)
    ed_arr =  LinRange(0.9*U, 0.75*U, 100) #u/2 +- u/2

    for i in eachindex(nsites_arr)
        Chem_Pot, Npart, E_gs, Em2_gs, Ep2_gs = Track_N(nsites_arr[i], ed_arr, 1, U)

        df = DataFrame(Onsite_Energy = ed_arr,
                    Chemical_Potential = Chem_Pot,
                    N_particles = Npart,
                    E1 = E_gs,
                    E1m2 = Em2_gs,
                    E1p2 = Ep2_gs)
        
        CSV.write(filenames[i], df)

    end
end

Nsites = [60, 100]#[10, 20, 40, 60, 100]
Filenames = ["U10N60_Track3.csv", "U10N100_Track3.csv"]#["U10N10_Track1.csv", "U10N20_Track1.csv", "U10N40_Track1.csv", "U10N60_Track1.csv", "U10N100_Track1.csv"]
u = 10

Roda_N(u, Nsites, Filenames)


#Nsites = [6]
#Filenames = ["PHSN6.csv"]
#u_arr = [0., 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 10., 12., 14., 16., 18., 20., 30., 40., 50.]
#u = 10

#ch = zeros(Float64, 0)
#Np = zeros(Float64, 0)
#for u in u_arr
#    ed_arr =  [u/2]#LinRange(1*u, 0*u, 100)
#
#    Chem_Pot1, Chem_Pot2, Npart = Track_N(Nsites[1], ed_arr, 1, u)
#    append!(ch, Chem_Pot1)
#    append!(Np, Npart)
    #print("\n", Chem_Pot1./u, "\n\n")
#end

#df = DataFrame(On_site_interaction = u_arr,
#                    Chemical_Potential = ch,
#                    N_eletrons = Np
#)
        
#CSV.write(Filenames[1], df)
#p = plot(u_arr, ch)
#display(p)
#p = scatter(Npart./Nsites[1], Chem_Pot1./u)
#p = scatter!(Npart./Nsites[1], Chem_Pot2./u)
#display(p)

#AA


#Roda_N(u, Nsites, Filenames)