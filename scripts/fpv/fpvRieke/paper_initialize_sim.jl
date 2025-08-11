
using DrWatson
@quickactivate "MonolithicFEMVLFS"

include(srcdir("fpv","paper_beam_joint_src.jl"))

using WaveSpec
using Plots
using .Constants
using .Jonswap
using DataFrames
using TimerOutputs
using TickTock 

using .beam

function constant_spectrum(start, finish,threshold) 
    ω_low = collect(start:0.02:threshold)
    ω_high = collect(threshold:0.1:finish)[2:end]
    ω = vcat(ω_low,ω_high)
    return ω
end

#Set the parameters
#h_b = [0.25,0.5,1.0]
#h_b = [1.0]
h_b = [0.20]

#constant wave spectrum 
ω = constant_spectrum(0.2,10,3)
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
    "numFloat" => [1,2,5,10,25,50]    
)


defaultParams = Dict(
    "mesh_size" => 1.0,
    "nz" => 12, #Remember, that we are using order 4 elements
    "mesh_rz" => 1.1,
    "emptyTankRef" => datadir("fpv_202508","empt",
        "Empt_depth=22.8_length_beam=100_mesh_size=1.0.jld2")
)

# ## Setting used in Rieke (2025) Energy paper
# nz = 12
# mesh_rz = 1.1

#creates all possible combinations
dicts = dict_list(allparams)
tick()

# ## Empty Tank
# function makesim(d::Dict)    
#     runDict = merge(d, defaultParams)

#     RAO_η, η_ϕ = beam.run_beam(runDict)

#     fulld = copy(runDict)
#     fulld["RAO_η"] = RAO_η
#     fulld["η_ϕ"] = η_ϕ
#     return fulld
# end


## Not Empty Tank
function makesim(d::Dict)    

    runDict = merge(d, defaultParams)
    
    RAO_η, η_ϕ, da_wavePrb, ηdof_scaled,
        EI, massPerArea = beam.run_beam(runDict)

    fulld = copy(runDict)
    fulld["RAO_η"] = RAO_η            
    fulld["η_ϕ"] = η_ϕ    
    fulld["EI"] = EI
    fulld["massPerArea"] = massPerArea    
    fulld["da_wavePrb"] = da_wavePrb
    fulld["ηdof"] = ηdof_scaled
    return fulld
end


for (i, d) in enumerate(dicts)
    tick()
    result = makesim(d) # is a dict conainting the input and output of simulation
    #safesave("data/sims2011/"*savename(d, "jld2"), result)
    safesave(datadir("fpv_202508", savename(d, "jld2")), result)
    tock()
end 

tock()



