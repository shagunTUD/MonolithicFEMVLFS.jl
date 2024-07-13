using DrWatson
using Parameters
using WaveSpec
using .Constants
@quickactivate "MonolithicFEMVLFS.jl"

# Here you may include files from the source directory
include(srcdir("viscoel3D","mem3D_modes_para_rad_unstr_iter.jl"))

lresDir::String = "data/sims_202403/run3D/mem3D_modes_res_free_2c"

# Warm-up run
mᵨ = [0.005] #970.0*0.005/1025 #mass per unit area of membrane / ρw
Tᵨ = [48.78] #1e7 stress x 0.005 thick #T/ρw
params = Memb3D_modes.Memb_params(
  resDir = lresDir,
  mesh_file = "models/mesh0/meshWorks.msh")
Memb3D_modes.main(params, mᵨ, Tᵨ)

# Production run
lparams = Memb3D_modes.Memb_params(
  mesh_file = "models/meshB2c/meshWorks.msh",
  resDir = lresDir,
  order = 1,
  order2 = 2,
  vtk_output = true,

  H0 = 10, #m #still-water depth

  nωₙ = 11,
  errLim = 1e-3,
  maxIter = 15,

  dbgmsg = false,

  # Domain These values do not matter for rad bnd condition
  Ld = 10,
  LΩ = 60,
  WΩ = 40,
  x₀ = -30,
  xdᵢₙ = -20,
  xdₒₜ = 20  
)
Memb3D_modes.main(lparams, mᵨ, Tᵨ)
