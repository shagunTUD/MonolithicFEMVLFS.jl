using DrWatson
using Parameters
using Gridap
using Printf
using WaveSpec
using .Constants
using MonolithicFEMVLFS.Resonator
@quickactivate "MonolithicFEMVLFS.jl"

# include(srcdir("lrmmModal","memModes_lrmm_free_iter.jl"))
include(srcdir("lrmmModal","memSysModes_lrmm_free_iter.jl"))


# Directory for results
resDir::String = "data/sims_202508/mem_modes_free/"


# mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
# tfac = [0.05, 0.10, 0.25, 0.50, 0.75, 0.8]

# mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
# tfac = [0.8]

mfac = [0.9] 
tfac = [0.1] 

ωr = collect(1.0:0.1:5.0)
resMᵨ = ones(length(ωr))
resKᵨ = [ iM*iωr*iωr for (iM, iωr) in zip(resMᵨ, ωr) ]

H0 = 10


# Common parameters for all runs
paramsBase = Membrane_modes.MembLR_params(

  vtk_output = false,

  H0 = H0, #m #still-water depth
  Lm = 2*H0, #m
  
  # mfac = imfac,
  # tfac = itfac,
  # rS_by_ρw = rS,

  nωₙ = 7, #number of natural frequencies to compute

  # Iterative solution for wet natural frequencies
  αRelax = 0.8,
  maxIter = 20,

  # initial guess for natural frequencies
  ωn_guess = zeros(ComplexF64, 1) .+ 0.5
)


# Loop over all combinations of mfac, tfac and resonator parameters
for (iresMᵨ, iresKᵨ) in zip(resMᵨ, resKᵨ)  

  rS = Resonator.Single( iresMᵨ, iresKᵨ, 0.0, Point(30.0, 0.0) )
  resonatorName = "res_m=" * @sprintf("%0.2f", iresMᵨ) *
    "_k=" * @sprintf("%0.2f", iresKᵨ)*"_mem"

  for imfac in mfac
    for itfac in tfac

      # Case directory
      caseName = "mem_modes_ten" * @sprintf("%0.2f", itfac) *
        "_mass" * @sprintf("%0.2f", imfac)

      caseDir = resDir*caseName
      # if( isdir(caseDir) )
      #   rm(caseDir, recursive=true) #remove old data
      #   @printf("Removed old data in %s\n", caseDir)
      #   # return
      # end
      # mkdir(caseDir)
      

      # Update paramsBase for each run
      params = Membrane_modes.MembLR_params(
        paramsBase;

        resDir = caseDir,
        fileName = resonatorName,

        mfac = imfac,
        tfac = itfac,
        rS_by_ρw = rS
      )

      # Run case
      Membrane_modes.run_case(params)
    end
  end
end