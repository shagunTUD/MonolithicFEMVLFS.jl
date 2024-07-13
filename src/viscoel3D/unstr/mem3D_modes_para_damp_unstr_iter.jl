module Memb3D_modes

using Revise
using Gridap
using Plots
using DrWatson
using WaveSpec
using .Constants
using LinearAlgebra
using TickTock
using DataFrames
using Printf
using Parameters
using GridapGmsh


function run_case(params, imᵨ = 0.9, iTᵨ = 0.1)
  
  function run_freq(ω)

    k = dispersionRelAng(H0, ω; msg=false)
    @show ω, k

    # Damping
    μ₀ = 2.5
    μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1]-x₀)/Ld))
    μ₁ₒᵤₜ(x) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒₜ)/Ld))
    μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
    μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
    ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)

  
    # Weak form: ω dependent
    k22(ϕ,w) = ∫( ∇(w)⋅∇(ϕ) )dΩ +
        ∫( -w * im * k * ϕ )dΓin + ∫( -w * im * k * ϕ )dΓot
  
    c23(κ,w) = ∫( im*ω*w*κ )dΓfs + 
      ∫( im*ω*w*κ - μ₂ᵢₙ*κ*w )dΓd1  +
      ∫( im*ω*w*κ - μ₂ₒᵤₜ*κ*w )dΓd2
  
    c32(ϕ,u) = ∫( -im*ω*u*ϕ )dΓfs +
      ∫( -im*ω*u*ϕ + μ₁ᵢₙ*∇ₙ(ϕ)*u )dΓd1 +  
      ∫( -im*ω*u*ϕ + μ₁ₒᵤₜ*∇ₙ(ϕ)*u )dΓd2
    
    # Global matrices: ω dependent
    K22 = get_matrix(AffineFEOperator( k22, l2, U_Ω, V_Ω ))
    C23 = get_matrix(AffineFEOperator( c23, l2, U_Γκ, V_Ω ))
    C32 = get_matrix(AffineFEOperator( c32, l3, U_Ω, V_Γκ ))
  
    # Solution
    tick()
    Mϕ = K22 - ( C23 * (Matrix(K33) \ C32) )
    Mhat = C12 * (Mϕ \ C21)
    Mtot = M11 + Mhat
    tock()
  
    # Eigen values
    # λ = LinearAlgebra.eigvals(Mtot\Matrix(K11))
    # V = LinearAlgebra.eigvecs(Mtot\Matrix(K11))  
    # # @show real.(λ[1:nωₙ])
    # # ωₙ = sqrt.(real.(λ))
    Ur, S, Vr = svd(Mtot\Matrix(K11))
    λ = reverse(S)
    V = reverse(Vr, dims=2)
    rλ = real.(λ[1:nωₙ])
    @show rλ
    return(rλ[1:nωₙ], V[:,1:nωₙ])
      
  end

  @unpack resDir, mesh_file, order, order2, vtk_output = params
  caseName = "ten" * @sprintf("%0.2f", iTᵨ) *"_mass" * @sprintf("%0.4f", imᵨ)
  name::String = resDir*"/mem_modes_"*caseName
  filename = name*"/mem"

  # if( isdir(name) )
  #   return
  # end
  # mkdir(name)

  if( ~isdir(name) )
    mkdir(name)
  end  

  @unpack H0 = params
  ρw = 1025 #kg/m3 water
  @show H0  #m #still-water depth

  # Membrane parameters
  @show g #defined in .Constants
  mᵨ = imᵨ  #mass per unit area of membrane / ρw
  Tᵨ = iTᵨ #T/ρw

  # Excitation wave parameters
  ω = 1.0

  # Domain
  @unpack Ld, LΩ, WΩ = params 
  @unpack x₀, xdᵢₙ, xdₒₜ = params
  @show LΩ, WΩ
  println()
  #@show Ld*k/2/π
  #@show cosh.(k*H0*0.5)./cosh.(k*H0)


  # Mesh  
  model = DiscreteModelFromFile(mesh_file)  


  # # Triangulations
  Ω = Interior(model) #same as Triangulation()
  Γ = Boundary(model, tags=["FreeSurface","Membrane"]) #same as BoundaryTriangulation()
  Γin = Boundary(model, tags="Inlet")
  Γot = Boundary(model, tags="Outlet")
  Γcyl = Boundary(model, tags="Cylinder")
  Γbot = Boundary(model, tags="Bottom")


  # Auxiliar functions
  function is_damping1(xs) # Check if an element is inside the damping zone 1
    n = length(xs)
    x = (1/n)*sum(xs)
    (x₀ <= x[1] <= xdᵢₙ ) * ( x[3] ≈ H0)
  end
  function is_damping2(xs) # Check if an element is inside the damping zone 2
    n = length(xs)
    x = (1/n)*sum(xs)
    (xdₒₜ <= x[1] ) * ( x[3] ≈ H0)
  end

  # Masking and Beam Triangulation
  Γκ = Boundary(model, tags = "FreeSurface")
  xΓκ = get_cell_coordinates(Γκ)
  Γd1_to_Γ_mask = lazy_map(is_damping1, xΓκ)
  Γd2_to_Γ_mask = lazy_map(is_damping2, xΓκ)
  Γd1 = Triangulation(Γκ, findall(Γd1_to_Γ_mask))
  Γd2 = Triangulation(Γκ, findall(Γd2_to_Γ_mask))  
  Γm = Boundary(model, tags = "Membrane")
  Γη = Γm  
  Γfs = Triangulation(Γκ, findall(!, Γd1_to_Γ_mask .| Γd2_to_Γ_mask))


  # # Construct the tag for membrane boundary
  # Λmb = Boundary(Γm)
  # xΛmb = get_cell_coordinates(Λmb)
  # xΛmb_n1 = findall(model.grid_topology.vertex_coordinates .== xΛmb[1])
  # xΛmb_n2 = findall(model.grid_topology.vertex_coordinates .== xΛmb[2])
  # new_entity = num_entities(labels_Ω) + 1
  # labels_Ω.d_to_dface_to_entity[1][xΛmb_n1[1]] = new_entity
  # labels_Ω.d_to_dface_to_entity[1][xΛmb_n2[1]] = new_entity
  # add_tag!(labels_Ω, "mem_bnd", [new_entity])


  writevtk(model, filename*"_model")
  if vtk_output == true
    writevtk(Ω,filename*"_O")
    writevtk(Γ,filename*"_G")
    writevtk(Γm,filename*"_Gm")  
    writevtk(Γfs,filename*"_Gfs")
    # writevtk(Λmb,filename*"_Lmb")  
  end


  # Measures
  degree = 2*order
  degree2 = 2*order2
  dΩ = Measure(Ω,degree)
  dΓm = Measure(Γm,degree2)
  dΓd1 = Measure(Γd1,degree2)
  dΓd2 = Measure(Γd2,degree2)
  dΓfs = Measure(Γfs,degree2)
  dΓin = Measure(Γin,degree)
  dΓot = Measure(Γot,degree)
  # dΛmb = Measure(Λmb,degree)


  # Normals
  # @show nΛmb = get_normal_vector(Λmb)


  # Dirichlet Fnc
  gη(x) = ComplexF64(0.0)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  reffe2 = ReferenceFE(lagrangian,Float64,order2)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γκ, reffe2, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  # V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
  #   vector_type=Vector{ComplexF64},
  #   dirichlet_tags=["mem_bnd"]) #diri
  V_Γη = TestFESpace(Γη, reffe2, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  # U_Γη = TrialFESpace(V_Γη, gη) #diri
  U_Γη = TrialFESpace(V_Γη)


  # Weak form: Constant matrices
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,0.0,1.0)  
  m11(η,v) = ∫( mᵨ*v*η )dΓm
  k11(η,v) = ∫( v*g*η + Tᵨ*∇(v)⋅∇(η) )dΓm #+  
              # ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri

  c12(ϕ,v) = ∫( v*ϕ )dΓm

  c21(η,w) = ∫( w*η )dΓm  

  k33(κ,u) = ∫( u*g*κ )dΓfs + 
    ∫(u*g*κ)dΓd1 + 
    ∫(u*g*κ)dΓd2

  l1(v) = ∫( 0*v )dΓm
  l2(w) = ∫( 0*w )dΩ
  l3(u) = ∫( 0*u )dΓfs 
  println("[MSG] Done Weak form")

  # Global matrices: constant matrices
  M11 = get_matrix(AffineFEOperator( m11, l1, U_Γη, V_Γη ))
  K11 = get_matrix(AffineFEOperator( k11, l1, U_Γη, V_Γη ))
  C12 = get_matrix(AffineFEOperator( c12, l1, U_Ω, V_Γη ))

  C21 = get_matrix(AffineFEOperator( c21, l2, U_Γη, V_Ω ))

  K33 = get_matrix(AffineFEOperator( k33, l3, U_Γκ, V_Γκ ))
  println("[MSG] Done Global matrices")
  println(K11 == transpose(K11))

  # #xp = range(xm₀, xm₁, size(V,2)+2)

  @unpack nωₙ, errLim, maxIter = params
  da_ωₙ = zeros(Float64, 1, nωₙ)
  @show ωₙ=zeros(Float64, 1, nωₙ) .+ ω
  da_V = []

  # # For index=1 not looping coz ωₙ[1] = 0.0
  # i = 1
  # ωₙ, V = run_freq(ω)
  # da_ωₙ[i] = ωₙ[1]
  # push!(da_V, V[:,i])

  for i in 1:nωₙ
    # global da_ωₙ, da_V  
    # global ωₙ, ω
    local V, lIter
    lIter = 0    
    Δω = 1
    ω = ωₙ[i]    
    while ((Δω > errLim) && (lIter < maxIter))
      # global ω, ωₙ
      rλ, V = run_freq(ω)
      ωₒ = ω      
      ωᵣ = sqrt(rλ[i])
      # if(i==1)
      #   #ω = 0.2 * ωₙ[i] + 0.8*ω
      #   ω = 0.0
      #   Δω = 0.0
      #   V = V*0.0
      # elseif(i==4)
      #   ω = 0.2 * ωᵣ + 0.8*ωₒ
      #   Δω = abs(ω - ωₒ)/ωₒ
      # else
      #   ω = 0.5 * ωᵣ + 0.5*ωₒ
      #   Δω = abs(ω - ωₒ)/ωₒ
      # end      
      ω = 0.5 * ωᵣ + 0.5*ωₒ
      Δω = abs(ω - ωₒ)/ωₒ
      # @show ωₙ      
      lIter += 1
      @show i, ω, Δω, lIter
    end
    da_ωₙ[i] = ω
    push!(da_V, V[:,i])
  end

  println(da_ωₙ)

  # xp = range(xm₀, xm₁, length(da_V[1]))

  data = Dict(
    "ωₙ" => da_ωₙ,
    "V" => da_V  
  )  

  wsave(filename*"_modesdata.jld2", data)
end


function main(params, mᵨ, Tᵨ)

  # mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
  # tfac = [0.05, 0.10, 0.25, 0.50, 0.75]


  # mᵨ = [0.005] #970.0*0.005/1025 #mass per unit area of membrane / ρw
  # Tᵨ = [48.78] #1e7 stress x 0.005 thick #T/ρw

  for imᵨ in mᵨ
    for iTᵨ in Tᵨ
      run_case(params, imᵨ, iTᵨ)
    end
  end
end


"""
Memb_params

Parameters for the VIV.jl module.
"""
@with_kw struct Memb_params

  mesh_file::String = "data/sims_202311/northHoyle_mesh/mesh/mesh0/meshWorks.msh"
  resDir::String = "data/sims_202309/mem3D_modes_res_free"
  order::Int = 1
  order2::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  nωₙ = 2
  errLim = 1e-1
  maxIter = 15

  dbgmsg = false

  # Domain 
  Ld = 20
  LΩ = 80
  WΩ = 40
  x₀ = -40
  xdᵢₙ = -20
  xdₒₜ = 20  

end


end