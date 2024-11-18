# nSpatialDims = 1
# Ngrid = 401
# dims  = fill(Ngrid, nSpatialDims+1)


# h₀s = collect(0.0001:0.0002:1.0)

# Ωperp = 10000    # Dimensional lumen footprint area
# Ωs    = h₀s.*Ωperp      # Dimensional lumen volume 


# h₀    = h₀s[40]
# Ω     = Ωs[40]      # Dimensional lumen volume 
# # Ωperp = Ω/h₀    # Dimensional lumen footprint area
# N     = 100     # Maximum polymer length 
# k_Cd  = 1.0 # Complex desorption rate
# k_Ca  = 0.01 # Complex adsorption rate
# k_Sd  = 1.0 # Substrate desorption rate
# k_Sa  = 0.2 # Substrate adsorption rate
# k₁    = 1.0   # Complex formation forward reaction rate 
# k₂    = 0.01   # Complex dissociation reverse reaction rate 
# k₃    = 0.1   # Product formation
# k₄    = 0.5  # Product dissociation 
# 𝓒     = 100000.0
# 𝓢     = 100000.0
# 𝓔     = 0.0001
# D_C   = 0.0000001  # Monomer/polymer diffusivity
# D_S   = 0.0000001  # Substrate diffusivity
# Tᵣstar= 10000000000.0  # Release time
# ϕ     = 0.5



nSpatialDims = 1
Ngrid = 401
dims  = fill(Ngrid, nSpatialDims+1)

# h₀s1 = collect(1e-5:2e-5:5e-4)
# h₀s2 = collect(5e-5+1e-6:1e-4:1e-2)
# h₀s = [h₀s1..., h₀s2...]

# h₀s2 = collect(1)

Ωperp = 10000    # Dimensional lumen footprint area
# Ωs    = h₀s.*Ωperp      # Dimensional lumen volume 

N     = 100     # Maximum polymer length 
k_Cd  = 1000000.0 # Complex desorption rate
k_Ca  = 0.01 # Complex adsorption rate
k_Sd  = 100000.0 # Substrate desorption rate
k_Sa  = 0.01 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.01   # Complex dissociation reverse reaction rate 
k₃    = 0.1   # Product formation
k₄    = 0.5  # Product dissociation 
𝓒     = 100000.0
𝓢     = 100000.0
𝓔     = 0.0001
D_C   = 0.0000001  # Monomer/polymer diffusivity
D_S   = 0.0000001  # Substrate diffusivity
Tᵣstar= 5000000000.0  # Release time
ϕ     = 0.5

λ = 0.05
σ = 0.2
