

function homogeneousWidthC(K₂, K₄, ttilde, β, α_C) 
    Etilde = K₂/(1+K₂)
    p1 = (1+α_C)/(π*Etilde*K₂*K₄*ttilde)
    p2 = ν*(1+α_C)-Etilde*β*ttilde
    p3 = 4*Etilde*K₂*K₄*(1+α_C)*ttilde
    return sqrt(p1)*exp(-p2^2/p3)
end



Ωperp = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd  = 200.0 # Complex desorption rate
k_Ca  = 2.0 # Complex adsorption rate
k_Sd  = 200.0 # Substrate desorption rate
k_Sa  = 1.1 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.1   # Complex dissociation reverse reaction rate 
k₃    = 1.0   # Product formation
k₄    = 1.0  # Product dissociation 
E_0   = 0.001
𝓒     = 100.0
𝓢     = 1000.0
D_C   = 0.01  # Monomer/polymer diffusivity
D_S   = 0.01  # Substrate diffusivity
Tᵣstar= 100.0  # Release time
ϕ     = 0.5

Nghost= 1           # Number of ghost points on each side of the domain 
Ngrid = 51

xMax = 100.0
xs   = collect(range(0.0, xMax, Ngrid+2*Nghost)) # Positions of discretised vertices in space

# h₀s = collect(0.1:0.1:3.0)
h₀s = collect(0.001:0.02:0.2001)