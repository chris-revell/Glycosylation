using Symbolics
using CairoMakie
using Latexify

@variables α_C α_S k_Ca k_Cd k_Sa k_Sd Ω 𝒜 h₀ L₀ k₁ k₂ k₃ k₄ K₂ K₃ K₄ C_b 𝒞 C_0 S_b 𝓢 S_0 σ E_0 𝓔 ϵ β N δ_C D_C δ_S D_S 𝓓 Tᵣ Tᵣstar Δ ϕ 𝓟star

lhs(e) = e.lhs
rhs(e) = e.rhs

expressions = [
    C_b ~ 𝒞/Ω,
    S_b ~ 𝓢/Ω,
    Ω ~ 𝒜*h₀,
    α_C ~ (k_Cd*Ω)/(2*k_Ca*𝒜),
    α_S ~ (k_Sd*Ω)/(2*k_Sa*𝒜),
    K₂ ~ (k₂/(k₁*C_b))*((2*k_Ca*𝒜 + k_Cd*Ω)/(k_Ca*Ω)),
    K₃ ~ k₃/k₁,
    K₄ ~ k₄/k₁,
    L₀ ~ sqrt(π)*Ω / (𝒜)^(1.5),
    C_0 ~ C_b*h₀/(2*(1+α_C)),
    S_0 ~ S_b*h₀/(2*(1+α_S)),
    σ ~ (k_Sa*S_b*(2*k_Ca*𝒜 + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*𝒜 + k_Sd*Ω)),
    𝓔 ~ 2*E_0*𝒜,
    σ ~ (k_Sa*S_b*(2*k_Ca*𝒜 + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*𝒜 + k_Sd*Ω)),
    ϵ ~ 𝓔*(2*k_Ca*𝒜 + k_Cd*Ω) / (2*k_Ca*C_b*𝒜),
    β ~ N*(σ*K₃ - K₂*K₄),
    δ_C ~ π*D_C/(k₁*𝓔),
    δ_S ~ π*D_S/(k₁*𝓔),
    𝓓 ~ α_C*δ_C*N^2*(K₂ + σ*K₃),
    Tᵣ ~ (k₁*𝓔*Tᵣstar)/(2*𝒜),
    β ~ N*(σ*K₃ - K₂*K₄),
    Δ ~ k₁*𝒞/(2*k₂*𝒜),
]

expressionsDict = Dict(lhs.(expressions).=>rhs.(expressions))

𝓟star = π/(2*ϕ) * (α_C*𝒞)/((1+α_C)^2) * (k₁*𝓔)/(2*𝒜) * K₂/(1+K₂) * (σ*K₃-K₂*K₄)/(N*(K₂+σ*K₃)) * (1/Tᵣ)
# expr = 𝓟star ~ π/(2*ϕ) * (α_C*𝒞)/((1+α_C)^2) * (k₁*𝓔)/(2*𝒜) * K₂/(1+K₂) * (σ*K₃-K₂*K₄)/(N*(K₂+σ*K₃)) * (1/Tᵣ)

sub1 = substitute(𝓟star, expressionsDict)
sub2 = substitute(sub1, expressionsDict)
sub3 = substitute(sub2, expressionsDict)

args = [h₀ 𝒜 N k₁ k₂ k₃ k₄ 𝓢 k_Sa k_Sd 𝒞 k_Ca k_Cd ϕ Tᵣstar]

𝓟starFunc = eval(build_function(sub3, args...))


valuesDict = Dict(
    𝒜 => 100.0,
    N     => 100,
    k_Cd => 0.9,
    k_Ca => 1.1,
    k_Sd => 0.9,
    k_Sa => 1.1,
    k₁   => 1.0,
    k₂   => 0.6,
    k₃   => 1.1,
    k₄   => 0.6,
    E_0 => 1.0,
    𝒞 => 100.0,
    𝓢 => 100.0,
    D_C  => 1.0,
    D_S  => 1.0,
    Tᵣstar  => 50.0,
    ϕ => 0.5,
)

𝒜 = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd = 0.9 # Complex desorption rate
k_Ca = 1.1 # Complex adsorption rate
k_Sd = 0.9 # Substrate desorption rate
k_Sa = 1.1 # Substrate adsorption rate
k₁   = 1.0   # Complex formation forward reaction rate 
k₂   = 0.6   # Complex dissociation reverse reaction rate 
k₃   = 1.1   # Product formation
k₄   = 0.6  # Product dissociation 
E_0 = 1.0
𝒞 = 100.0
𝓢 = 100.0
D_C  = 1.0  # Monomer/polymer diffusivity
D_S  = 1.0  # Substrate diffusivity
Tᵣstar  = 100.0  # Release time
ϕ = 0.5


𝓟starFunc(1.1, 𝒜, N, k₁, k₂, k₃, k₄, 𝓢, k_Sa, k_Sd, 𝒞, k_Ca, k_Cd, ϕ, Tᵣstar)


# 𝓟starVal = substitute(sub3, valuesDict)

h₀s = collect(0.1:0.1:3.0)

# Ps = zeros(length(h₀s))
Ps = 𝓟starFunc.(h₀s, 𝒜, N, k₁, k₂, k₃, k₄, 𝓢, k_Sa, k_Sd, 𝒞, k_Ca, k_Cd, ϕ, Tᵣstar)


for i=1:length(h₀s)
    num = substitute(𝓟starVal, h₀=>h₀s[i])
    Ps[i] = num.val
end

fig = CairoMakie.Figure(size=(500,500))
ax = Axis(fig[1,1])
ax.xlabel = "h₀"
ax.ylabel = L"𝓟^*"
ylims!(ax, (0.0,maximum(Ps)))
xlims!(ax, (0.0,maximum(h₀s)))
lines!(ax, h₀s, Ps)
display(fig)

save("parameterSweep.png",fig)


#%%

expressions2 = []
for i=1:length(expressions)
    innr = [expressions[i].rhs]
    for exp in expressions
        innr[1] = substitute(innr[1], Dict(lhs.(expressions).=>rhs.(expressions)))
    end
    push!(expressions2, expressions[i].lhs~innr[1])
end

# substitute.(expressions2, valuesDict)
tst = [substitute(exp, valuesDict) for exp in expressions2]
#%%


render(latexify(sub3))



#%%

# expr1 = Ω ~ 𝒜*h₀
# expr2 = α_C ~ (k_Cd*Ω)/(2*k_Ca*𝒜)                          # Balance of complex in bulk to complex on membrane
# expr3 = α_S ~ (k_Sd*Ω)/(2*k_Sa*𝒜)                          # Balance of substrate in bulk to substrate on membrane
# expr4 = K₂ ~ (k₂/(k₁*C_b))*((2*k_Ca*𝒜 + k_Cd*Ω)/(k_Ca*Ω)) # Non-dimensionalised complex formation net reaction rate k₂/(k₁*C_0)
# expr5 = K₃ ~ k₃/k₁                                            # Non-dimensionalised product formation rate
# expr6 = K₄ ~ k₄/k₁                                            # Non-dimensionalised prodict dissociation rate
# # expr h₀ ~ Ω/𝒜                                          # Mean thickness 
# expr7 = L₀ ~ sqrt(π)*Ω / (𝒜)^(1.5)                        # Mean radius 
# # expr C_b~ 𝒞/Ω                     # Initial monomer bulk concentration 
# # expr S_b~ 𝓢/Ω                     # Initial substrate mass
# expr8 = C_0 ~ C_b*h₀/(2*(1+α_C))        # Early surface monomer concentration
# expr9 = S_0 ~ S_b*h₀/(2*(1+α_S))      # Early surface substrate concentration 
# # expr E_0~ 𝓔/(2*𝒜)             # Total enzyme mass
# expr10 = 𝒞  ~ C_b*Ω                     # Initial monomer mass
# expr11 = σ  ~ (k_Sa*S_b*(2*k_Ca*𝒜 + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*𝒜 + k_Sd*Ω))
# expr12 = 𝓢  ~ S_b*Ω                     # Initial substrate mass
# expr13 = 𝓔  ~ 2*E_0*𝒜               # Total enzyme mass
# expr14 = σ  ~ (k_Sa*S_b*(2*k_Ca*𝒜 + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*𝒜 + k_Sd*Ω)) # S_0/C_0
# expr15 = ϵ  ~ 𝓔*(2*k_Ca*𝒜 + k_Cd*Ω) / (2*k_Ca*C_b*𝒜) # E_0/C_0
# expr16 = β  ~ N*(σ*K₃ - K₂*K₄)
# expr17 = δ_C ~ π*D_C/(k₁*𝓔)
# expr18 = δ_S ~ π*D_S/(k₁*𝓔)
# expr19 = 𝓓  ~ α_C*δ_C*N^2*(K₂ + σ*K₃)
# expr20 = Tᵣ ~ (k₁*𝓔*Tᵣstar)/(2*𝒜)
# expr21 = β ~ N*(σ*K₃ - K₂*K₄)
# expr22 = Δ  ~ k₁*𝒞/(2*k₂*𝒜)

# Ω ~ 𝒜*h₀
# α_C ~ (k_Cd*Ω)/(2*k_Ca*𝒜)                          # Balance of complex in bulk to complex on membrane
# α_S ~ (k_Sd*Ω)/(2*k_Sa*𝒜)                          # Balance of substrate in bulk to substrate on membrane
# K₂ ~ (k₂/(k₁*C_b))*((2*k_Ca*𝒜 + k_Cd*Ω)/(k_Ca*Ω)) # Non-dimensionalised complex formation net reaction rate k₂/(k₁*C_0)
# K₃ ~ k₃/k₁                                            # Non-dimensionalised product formation rate
# K₄ ~ k₄/k₁                                            # Non-dimensionalised prodict dissociation rate
# # expr h₀ ~ Ω/𝒜                                          # Mean thickness 
# L₀ ~ sqrt(π)*Ω / (𝒜)^(1.5)                        # Mean radius 
# # expr C_b~ 𝒞/Ω                     # Initial monomer bulk concentration 
# # expr S_b~ 𝓢/Ω                     # Initial substrate mass
# C_0 ~ C_b*h₀/(2*(1+α_C))        # Early surface monomer concentration
# S_0 ~ S_b*h₀/(2*(1+α_S))      # Early surface substrate concentration 
# # expr E_0~ 𝓔/(2*𝒜)             # Total enzyme mass
# 𝒞  ~ C_b*Ω                     # Initial monomer mass
# σ  ~ (k_Sa*S_b*(2*k_Ca*𝒜 + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*𝒜 + k_Sd*Ω))
# 𝓢  ~ S_b*Ω                     # Initial substrate mass
# 𝓔  ~ 2*E_0*𝒜               # Total enzyme mass
# σ  ~ (k_Sa*S_b*(2*k_Ca*𝒜 + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*𝒜 + k_Sd*Ω)) # S_0/C_0
# ϵ  ~ 𝓔*(2*k_Ca*𝒜 + k_Cd*Ω) / (2*k_Ca*C_b*𝒜) # E_0/C_0
# β  ~ N*(σ*K₃ - K₂*K₄)
# δ_C ~ π*D_C/(k₁*𝓔)
# δ_S ~ π*D_S/(k₁*𝓔)
# 𝓓  ~ α_C*δ_C*N^2*(K₂ + σ*K₃)
# Tᵣ ~ (k₁*𝓔*Tᵣstar)/(2*𝒜)
# β ~ N*(σ*K₃ - K₂*K₄)
# Δ  ~ k₁*𝒞/(2*k₂*𝒜)


# @variables α_C α_S k_Ca k_Cd k_Sa k_Sd Ω 𝒜 h₀ L₀ k₁ k₂ k₃ k₄ K₂ K₃ K₄ C_b 𝒞 C_0 S_b 𝓢 S_0 σ E_0 𝓔 ϵ β N δ_C D_C δ_S D_S 𝓓 Tᵣ Tᵣstar Δ ϕ 𝓟star

# Ω = 𝒜*h₀
# C_b = 𝒞/Ω
# S_b = 𝓢/Ω 
# α_C = (k_Cd*Ω)/(2*k_Ca*𝒜)                          
# α_S = (k_Sd*Ω)/(2*k_Sa*𝒜)                          
# K₂ = (k₂/(k₁*C_b))*((2*k_Ca*𝒜 + k_Cd*Ω)/(k_Ca*Ω)) 
# K₃ = k₃/k₁                                            
# K₄ = k₄/k₁                                            
# # expr h₀= Ω/𝒜                                    
# L₀ = sqrt(π)*Ω / (𝒜)^(1.5)                        
# # expr C_b~ 𝒞/Ω                     
# # expr S_b~ 𝓢/Ω                     
# C_0 = C_b*h₀/(2*(1+α_C))        
# S_0 = S_b*h₀/(2*(1+α_S))      
# # expr E_0~ 𝓔/(2*𝒜)             
# # 𝒞 = C_b*Ω                    
# σ = (k_Sa*S_b*(2*k_Ca*𝒜 + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*𝒜 + k_Sd*Ω))
# # 𝓢 = S_b*Ω                     
# 𝓔 = 2*E_0*𝒜               
# σ = (k_Sa*S_b*(2*k_Ca*𝒜 + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*𝒜 + k_Sd*Ω)) # S_0/C_0
# ϵ = 𝓔*(2*k_Ca*𝒜 + k_Cd*Ω) / (2*k_Ca*C_b*𝒜) # E_0/C_0
# β = N*(σ*K₃ - K₂*K₄)
# δ_C = π*D_C/(k₁*𝓔)
# δ_S = π*D_S/(k₁*𝓔)
# 𝓓 = α_C*δ_C*N^2*(K₂ + σ*K₃)
# Tᵣ = (k₁*𝓔*Tᵣstar)/(2*𝒜)
# β = N*(σ*K₃ - K₂*K₄)
# Δ = k₁*𝒞/(2*k₂*𝒜)


# 𝓟starSimpl = simplify(sub1)


# expressions0 = [
#     K₃ ~ k₃/k₁,
#     K₄ ~ k₄/k₁,
#     𝓔 ~ 2*E_0*𝒜,
#     Δ ~ k₁*𝒞/(2*k₂*𝒜),
#     Ω ~ 𝒜*h₀,
# ]

# expressions1 = [
#     α_C ~ (k_Cd*Ω)/(2*k_Ca*𝒜),
#     α_S ~ (k_Sd*Ω)/(2*k_Sa*𝒜),
#     L₀ ~ sqrt(π)*Ω / (𝒜)^(1.5),
#     δ_C ~ π*D_C/(k₁*𝓔),
#     δ_S ~ π*D_S/(k₁*𝓔),
#     Tᵣ ~ (k₁*𝓔*Tᵣstar)/(2*𝒜),
#     C_b ~ 𝒞/Ω,
#     S_b ~ 𝓢/Ω,
# ]

# expressions2 = [
#     K₂ ~ (k₂/(k₁*C_b))*((2*k_Ca*𝒜 + k_Cd*Ω)/(k_Ca*Ω)),
#     C_0 ~ C_b*h₀/(2*(1+α_C)),
#     S_0 ~ S_b*h₀/(2*(1+α_S)),
#     σ ~ (k_Sa*S_b*(2*k_Ca*𝒜 + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*𝒜 + k_Sd*Ω)),
#     ϵ ~ 𝓔*(2*k_Ca*𝒜 + k_Cd*Ω) / (2*k_Ca*C_b*𝒜), # E_0/C_0
#     β ~ N*(σ*K₃ - K₂*K₄),
#     𝓓 ~ α_C*δ_C*N^2*(K₂ + σ*K₃),
# ]
# @variables k_Ca k_Cd k_Sa k_Sd 𝒜 h₀ k₁ k₂ k₃ k₄ C_b S_b E_0 N D_C D_S Tᵣstar ϕ
# @variables Ω(𝒜,h₀) α_C(k_Cd,Ω,k_Ca,𝒜) α_S(Ω,k_Sa,𝒜) L₀(Ω,𝒜) 
# @variables K₂(k₁,k₂,C_b,k_Ca,𝒜,k_Cd,Ω) K₃(k₁,k₃) K₄(k₁,k₄) 𝒞(C_b,Ω) C_0(C_b,h₀,α_C) 𝓢(S_b,Ω) S_0(S_b,h₀,α_S) 
# @variables σ(S_b,k_Cd,k_Ca,C_b,k_Sa,k_Sd,𝒜,Ω) 𝓔(E_0,𝒜) 
# @variables ϵ(𝓔,k_Cd,k_Ca,C_b,Ω,𝒜) β(N,σ,K₃,K₂,K₄) 
# @variables δ_C(D_C,k₁,𝓔) δ_S(D_S,k₁,𝓔) 𝓓(α_C,δ_C,N,K₂,σ,K₃) Tᵣ(k₁,𝓔,Tᵣstar,𝒜) Δ(k₁,𝒞,k₂,𝒜) 
# @variables 𝓟star(ϕ,𝒞,α_C,k₁,𝓔,𝒜,K₂,K₄,N,K₂,σ,K₃,Tᵣ)
