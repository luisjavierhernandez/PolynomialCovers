using PolynomialCovers

# Polynomial P(z) = z^2 - 1
P = PolynomialCover([1.0,0.0,-1.0])

# Basepoint in the codomain
b = basepoint(P)

# Compute the fiber
F = fiber(P)

println("Initial fiber:")
println(F)

# Construct a loop around the origin (critical value)
γ = basepoint_initial_loop(P)

# Transport the fiber along the loop
transport!(F, γ)

println("Fiber after lifting the loop:")
println(F)