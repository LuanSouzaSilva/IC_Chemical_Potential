using Pkg

#=
Pkg.add("ITensors")
Pkg.add("Plots")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("SparseArrays")
Pkg.add("Arpack")
Pkg.add("LinearAlgebra")
Pkg.add("ProgressBars")
=#

include("DMRG_Hubbard.jl")

#Pkg.instantiate()

using ITensors, Plots, SparseArrays, Arpack, CSV, DataFrames, LinearAlgebra, ProgressBars

JULIA_NUM_THREADS=12

Threads.nthreads()

ITensors.enable_threaded_blocksparse()
NDTensors.Strided.disable_threads()
BLAS.set_num_threads(1)


nsites = 20

t = 1
U = 10

ed_arr =  LinRange(-1*u, 0, 50) #-u/2 +- u/2

GS_trial, Npart_trial_exp = Hubbard_DMRG(nsites, t, U, ed_arr[1])
Npart_trial = round(sum(Npart_trial_exp))

print(GS_trial, "      ", Npart_trial, "\n\n")

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

Npart_arr = zeros(Float64, 0)

append!(Npart_arr, Npart)
for  i in eachindex(ed_arr)[2:length(ed_arr)]
    local GS_Nconserve = Npart_DMRG(nsites, t, U, ed_arr[i], Npart)
    local GS_Nconservep1 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart + 1)
    local GS_Nconservem1 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart - 1)


    if GS_Nconserve < GS_Nconservem1 && GS_Nconserve < GS_Nconservep1
        global Npart += 0
    elseif GS_Nconservem1 < GS_Nconserve
        global Npart -= 1
    elseif GS_Nconservep1 < GS_Nconserve
        global Npart += 1
    end
    
    print(i, " -- ", ed_arr[i], " -- ", Npart, "\n")
    #print(GS_Nconserve, "        ", GS_Nconservem1, "         ", GS_Nconservep1, "\n")

    append!(Npart_arr, Npart)

end


Chem_pot = zeros(Float64, 0)
for i in eachindex(ed_arr)
    GSp1 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart_arr[i]+1)
    GSm1 = Npart_DMRG(nsites, t, U, ed_arr[i], Npart_arr[i]-1)

    mu = (GSp1 - GSm1)/2

    append!(Chem_pot, mu)
end

p = scatter(1/2 .+ ed_arr./U, Chem_pot./U)
display(p)

#AA