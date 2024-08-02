using Catalyst

@parameters kSa kSd kCa kCd k₁ k₂ k₃ k₄
t = default_t()
@species S(t) Sₛ(t) Cₙ(t) Cₙₛ(t) Cₙ₊₁(t) Cₙ₊₁ₛ(t) E(t) Qₙ(t)

reactionsVec = Reaction[]

push!(reactionsVec, Reaction(kSa, [S], [Sₛ], [1], [1]))
push!(reactionsVec, Reaction(kSd, [Sₛ], [S], [1], [1]))
push!(reactionsVec, Reaction(kCa, [Cₙ], [Cₙₛ], [1], [1]))
push!(reactionsVec, Reaction(kCd, [Cₙₛ], [Cₙ], [1], [1]))
push!(reactionsVec, Reaction(k₁, [Cₙₛ, E], [Qₙ], [1, 1], [1]))
push!(reactionsVec, Reaction(k₂, [Qₙ], [Cₙₛ, E], [1], [1, 1]))
push!(reactionsVec, Reaction(k₃, [Qₙ, Sₛ], [Cₙ₊₁ₛ, E], [1, 1], [1, 1]))
push!(reactionsVec, Reaction(k₄, [Cₙ₊₁ₛ, E], [Qₙ, Sₛ], [1, 1], [1, 1]))
push!(reactionsVec, Reaction(kSd, [Cₙ₊₁ₛ], [Cₙ₊₁], [1, 1], [1, 1]))
push!(reactionsVec, Reaction(kSa, [Cₙ₊₁], [Cₙ₊₁ₛ], [1, 1], [1, 1]))


# rx1 = Reaction(kSa, [S], [Sₛ], [1], [1])
# rx2 = Reaction(kSd, [Sₛ], [S], [1], [1])
# rx3 = Reaction(kCa, [Cₙ], [Cₙₛ], [1], [1])
# rx4 = Reaction(kCd, [Cₙₛ], [Cₙ], [1], [1])

# rx5 = Reaction(k₁, [Cₙₛ, E], [Qₙ], [1, 1], [1])
# rx6 = Reaction(k₂, [Qₙ], [Cₙₛ, E], [1], [1, 1])
# rx7 = Reaction(k₃, [Qₙ, Sₛ], [Cₙ₊₁ₛ, E], [1, 1], [1, 1])
# rx8 = Reaction(k₄, [Cₙ₊₁ₛ, E], [Qₙ, Sₛ], [1, 1], [1, 1])

@named network = ReactionSystem(reactionsVec, t)
networkComplete = complete(network)
Graph(networkComplete)
complexgraph(networkComplete)