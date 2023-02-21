#module Membrane_modes

using Revise
using Gridap
using Plots
using DrWatson
using WaveSpec
using .Constants
using LinearAlgebra


name::String = "data/sims_202302/mem_modes"
order::Int = 1
vtk_output::Bool = true
filename = name*"/mem"

ρw = 1025 #kg/m3 water
H0 = 10 #m #still-water depth

# Membrane parameters
@show Lm = 2*H0 #m
@show g #defined in .Constants
@show mᵨ = 0.9 #mass per unit area of membrane / ρw
@show Tᵨ = 0.1*g*H0*H0 #T/ρw
@show τ = 0.05#damping coeff


# Domain 
nx = 200
ny = 20
mesh_ry = 1.1 #Ratio for Geometric progression of eleSize
LΩ = 2*H0
x₀ = 0
domain =  (x₀, x₀+LΩ, -H0, 0.0)
partition = (nx, ny)
@show LΩ
@show domain
@show partition
@show LΩ/nx
@show H0/ny
println()


# Mesh
function f_y(y, r, n, H0)
  # Mesh along depth as a GP
  # Depth is 0 to -H0
  if(r ≈ 1.0)
    return y  
  else
    a0 = H0 * (r-1) / (r^n - 1)    
    if y ≈ 0
      return 0.0
    end
    j = abs(y) / H0 * n  
    return -a0 / (r-1) * (r^j - 1)
  end
end
map(x) = VectorValue( x[1], f_y(x[2], mesh_ry, ny, H0) )
model = CartesianDiscreteModel(domain,partition,map=map)


# Labelling
labels_Ω = get_face_labeling(model)
add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)
add_tag_from_tags!(labels_Ω, "water", [9])       # assign the label "water" to the entity 9 (interior)


# Triangulations
Ω = Interior(model) #same as Triangulation()
Γm = Boundary(model,tags="surface") #same as BoundaryTriangulation()
Γin = Boundary(model,tags="inlet")
Γot = Boundary(model,tags="outlet")


# Normals
nΓin = get_normal_vector(Γin)
nΓot = get_normal_vector(Γot)


writevtk(model, filename*"_model")
if vtk_output == true
  writevtk(Ω,filename*"_O")  
  writevtk(Γm,filename*"_Gm")    
end


# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓm = Measure(Γm,degree)
dΓin = Measure(Γin,degree)
dΓot = Measure(Γot,degree)


# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω, reffe, conformity=:H1,
  vector_type=Vector{ComplexF64})
V_Γη = TestFESpace(Γm, reffe, conformity=:H1,
  vector_type=Vector{ComplexF64})
U_Ω = TrialFESpace(V_Ω)
U_Γη = TrialFESpace(V_Γη)

ω = 2.4
k = dispersionRelAng(H0, ω)
@show tanh(k*H0)
Bc1 = 1/9.81/tanh(k*H0)
Bc2 = 1/9.81/tanh(k*H0)
# Weak form
mbb(η,v) = ∫( mᵨ*v*η )dΓm
kbb(η,v) = ∫( v*g*η + Tᵨ*∇(v)⋅∇(η) )dΓm
lb(v) = ∫( 0*v )dΓm

cwb(ϕ,v) = ∫( v*ϕ )dΓm

cbw(η,w) = ∫( w*η )dΓm

kww(ϕ,w) =  ∫(  ∇(w)⋅∇(ϕ) )dΩ - 
            im*ω^2*( ∫(Bc1*w*(∇(ϕ)⋅VectorValue(1.0,0.0)))dΓin 
              - ∫(Bc2*w*(∇(ϕ)⋅VectorValue(1.0,0.0)))dΓot )

lw(w) = ∫( 0*w )dΩ
# l((w,u,v)) =  ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1

Mbb = get_matrix(AffineFEOperator(mbb,lb,U_Γη,V_Γη))
Kbb = get_matrix(AffineFEOperator(kbb,lb,U_Γη,V_Γη))
Cwb = get_matrix(AffineFEOperator(cwb,lb,U_Ω,V_Γη))

Cbw = get_matrix(AffineFEOperator(cbw,lw,U_Γη,V_Ω))
Kww = get_matrix(AffineFEOperator(kww,lw,U_Ω,V_Ω))

Kww_inv = inv(Matrix(Kww))
Mhat = Cwb*Kww_inv*Cbw

λ = LinearAlgebra.eigvals(Matrix((Mbb+Mhat)\Kbb))
V = LinearAlgebra.eigvecs(Matrix((Mbb+Mhat)\Kbb))
@show sum(imag.(λ))
ω = sqrt.(abs.(λ))
@show ind = findall(ω.<5)
@show ω[ind]

xp = range(x₀, x₀+LΩ, nx*order+1)
for i in ind
  ωi = ω[i]
  Vi = V[:,i]
  plt1 = plot(xp, abs.(Vi), linewidth=3, 
    xlabel = "x (m)",
    ylabel = "Amp (m)", legend=false,
    title = "Mode_$i ω = $ωi")
  savefig(plt1,name*"/mode_$i.png")
end

#end