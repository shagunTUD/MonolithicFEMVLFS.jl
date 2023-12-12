
using DrWatson
@quickactivate "Offshore PV"

export give_material

function give_material(name::String)
    Materials_dict = Dict( "steel" => [932.26032, 2.15e11, 0.2192 ], #store material properties: density (kg/m^3) and young's modulus N/m^2 and in case of hollow material ration Ix_hollow/Ix_solid
                        "hdpe" => [250, 12e9,1],
                        "neopren" => [116, 560e3,1],
                        "cfrp" => [960.432, 40e9, 0.8704], #Carbon fiber reinforced polymer based on https://doi.org/10.1016/j.matpr.2020.07.636, lower end is used (conservative)
                        "gfrp" => [937.8, 35e9, 0.9375],  #Glass fiber
                        "eps" => [15, 6.5e6, 1],    #
                        "pvc" => [896.432, 52e6, 0.8704]
    )

    ρ = Materials_dict[name][1]
    E = Materials_dict[name][2]
    γ = Materials_dict[name][3]
    
    return ρ,E,γ
    #return keys(Materials)
end




