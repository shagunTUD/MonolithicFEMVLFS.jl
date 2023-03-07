## Natural frequency ω₀
mbb((η,θ₁,θ₂),(v,z₁,z₂ )) = ∫( d₀*η*v + (I/ρ_w)*(θ₁⋅z₁) )dΓb1 +
∫( d₀*η*v + (I/ρ_w)*(θ₂⋅z₂) )dΓb2
cbb((η,θ₁,θ₂),(v,z₁,z₂ )) = ∫( (jump_Λj(z₁,z₂) ⋅ nΛj1D) * (cᵣ/ρ_w * (jump_Λj(θ₁,θ₂) ⋅ nΛj1D)) )dΛj
kbb((η,θ₁,θ₂),(v,z₁,z₂ )) = ∫( ( ₓ∘(∇(v)) - z₁ ) ⋅ ( a₂[1]*(ₓ∘(∇(η))-θ₁) ) + ₓ∘(∇(z₁))⊙(a₁[1]*(ₓ∘(∇(θ₁))))+g*η*v)dΓb1 +
∫( ( ₓ∘(∇(v)) - z₂ ) ⋅ ( a₂[2] *(ₓ∘(∇(η))-θ₂) ) + ₓ∘(∇(z₂))⊙(a₁[2]*(ₓ∘(∇(θ₂))))+g*η*v)dΓb2 +
∫( (jump_Λj(z₁,z₂) ⋅ nΛj1D) * (kᵣ/ρ_w * (jump_Λj(θ₁,θ₂) ⋅ nΛj1D)) )dΛj
lb((v,z₁,z₂ )) = 0
cwb((ϕ,κ),(v,z₁,z₂)) = ∫( ϕ*v )dΓb1 +
∫( ϕ*v )dΓb2
cbw((η,θ₁,θ₂),(w,u)) = ∫( η*w )dΓb1 +
∫( η*w )dΓb2
kww((ϕ,κ),(w,u)) = ∫( ∇(w)⋅∇(ϕ) )dΩ +
∫( βₕ*u*g*κ )dΓfs +
∫( βₕ*u*g*κ)dΓd1 + #- μ₁ᵢₙ*κ*w - μ₂ᵢₙ*ϕ*w/g )dΓd1 +
∫( βₕ*u*g*κ)dΓd2 # - μ₁ₒᵤₜ*κ*w - μ₂ₒᵤₜ*ϕ*w/g )dΓd2
# cww((ϕ,κ),(w,u)) = ∫( -1im*ω*βₕ*ϕ*u - (-1im*ω)*κ*w + βₕ*αₕ*w*g*κ )dΓfs +
# ∫( -1im*ω*βₕ*ϕ*u - (-1im*ω)*κ*w + βₕ*αₕ*w*g*κ )dΓd1 +
# ∫( -1im*ω*βₕ*ϕ*u - (-1im*ω)*κ*w + βₕ*αₕ*w*g*κ )dΓd2
lw((w,u)) = 0
Mbb = get_matrix(AffineFEOperator(mbb,lb,Xb,Yb))
Cbb = get_matrix(AffineFEOperator(cbb,lb,Xb,Yb))
Kbb = get_matrix(AffineFEOperator(kbb,lb,Xb,Yb))
Cbw = get_matrix(AffineFEOperator(cbw,lw,Xb,Yw))
Cwb = get_matrix(AffineFEOperator(cwb,lb,Xw,Yb))
Kww = get_matrix(AffineFEOperator(kww,lw,Xw,Yw))
Kww_inv = inv(Matrix(Kww))
M_hat = Cwb*Kww_inv*-Cbw
@show sum(imag.(Mbb-M_hat))
@show sum(imag.(Kbb))
n = size(Mbb)[1]
#λ,V = eigs(Kbb,Mbb-M_hat,nev=(n))
λ = LinearAlgebra.eigvals(Matrix((Mbb-M_hat)\Kbb))
V = LinearAlgebra.eigvecs(Matrix((Mbb-M_hat)\Kbb))
@show sum(imag.(λ))
ω = sqrt.(abs.(λ))
f = ω/2π
