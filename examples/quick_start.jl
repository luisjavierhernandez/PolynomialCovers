using PolynomialCovers

# Define a polynomial cover
P = PolynomialCover([1, 0, -1])     # P(z) = z^2 - 1

# Compute a safe basepoint in the codomain
b = basepoint(P)

# Compute the fiber P⁻¹(b)
F = fiber(P)

println("Fiber points:")
println(F)

# Construct a path avoiding critical values
γ = safe_path_arc(b, -b)

# Transport the fiber along the path
transport!(F, γ)

println("Fiber after transport:")
println(F)