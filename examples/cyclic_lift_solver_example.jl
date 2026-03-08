using PolynomialCovers

P = PolynomialCover(1, 0, -1)

# Radius containing the critical values
S = critical_values_radius(P)

println("Critical values lie in the disk |z| ≤ $S")

# Choose a value far from the critical values
b = 3S

# Solve P(z) = b
z = cyclic_lift_solver(P, b)

println("Solution:", z)