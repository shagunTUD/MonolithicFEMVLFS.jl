
using DrWatson
@quickactivate "MonolithicFEMVLFS"

include(srcdir("fpv","paper_beam_joint_src_singleFreq.jl"))

using WaveSpec
using Plots
using .Constants
using .Jonswap
using DataFrames
using TimerOutputs
using TickTock 

using .beam


#Set the parameters
#h_b = [0.25,0.5,1.0]
#h_b = [1.0]
h_b = [0.20]

#constant wave spectrum 
ω = [1.5]
# ω = collect(0.2:1:3.2)
α = vec(repeat([0], length(ω),1))
A_w = vec(repeat([1], length(ω),1))


allparams = Dict(    
    "length_beam" => [100],
    "h_b" => h_b,                            
    # "material" => ["eps", "cfrp", "gfrp"],# ["eps", "cfrp", "gfrp", "neopren", "pvc", "hdpe", "steel"], #store youngs modulus and density 
    "material" => ["hdpe"],
    "phase"=> [α],
    "omega" => [ω],
    "depth" => [22.8],
    "amplitude" => [A_w],
    "numFloat" => [10],

    "mesh_size" => [1.0],
    "nz" => [12], #Remember, that we are using order 4 elements
    "mesh_rz" => [1.1]
)

# ## Setting used in Rieke (2025) Energy paper
# nz = 12
# mesh_rz = 1.1

#creates all possible combinations
dicts = dict_list(allparams)
tick()

## Empty Tank
# function makesim(d::Dict)    
#     RAO_η, RAO_ηx, η_ϕ, ηx_ϕ = beam.run_beam(d)
#     fulld = copy(d)
#     fulld["RAO_η"] = RAO_η
#     fulld["RAO_ηx"] = RAO_ηx
#     fulld["η_ϕ"] = η_ϕ
#     fulld["ηx_ϕ"] = ηx_ϕ
#     return fulld
# end

# ## Not Empty Tank
# function makesim(d::Dict)    
#     RAO_η, RAO_ηx, η_ϕ, ηx_ϕ, 
#         RAO_ηxx, RAO_ηxxx, da_wavePrb, ηdof_scaled,
#         EI, massPerArea = beam.run_beam(d)
#     fulld = copy(d)
#     fulld["RAO_η"] = RAO_η
#     fulld["RAO_ηx"] = RAO_ηx
#     fulld["RAO_ηxx"] = RAO_ηxx
#     fulld["RAO_ηxxx"] = RAO_ηxxx
#     fulld["η_ϕ"] = η_ϕ
#     fulld["ηx_ϕ"] = ηx_ϕ
#     fulld["EI"] = EI
#     fulld["massPerArea"] = massPerArea    
#     fulld["da_wavePrb"] = da_wavePrb
#     fulld["ηdof"] = ηdof_scaled
#     return fulld
# end

## Not Empty Tank
function makesim(d::Dict)    
    beam.run_beam(d, vtk_output=true, 
        filename = datadir("fpv_202508","run","beam"))    
end


for (i, d) in enumerate(dicts)
    tick()
    makesim(d) # is a dict conainting the input and output of simulation
    tock()
end 

tock()



