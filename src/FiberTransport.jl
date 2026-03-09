
# Estructura principal para la fibra
mutable struct Fiber{T<:Real}
	P::PolynomialCover{T}
	xs::Vector{Complex{T}}
	w::Complex{T}
	S::Real
end

# Construcción de la fibra usando monodromía en el infinito

"""
    fiber(P::PolynomialCover{T}; nl::Int=8, s::Real=3, tol::Real=1e-9) -> Fiber

Compute the fiber of a polynomial cover `P` at a generic basepoint using
path lifting along a loop around the "point at infinity" (monodromy).

# Arguments
- `P::PolynomialCover{T}`: The polynomial cover for which to compute the fiber.
- `nl::Int=8`: Number of linear segments in the basepoint loop.
- `s::Real=3`: Scaling factor for the initial basepoint in the codomain.
- `tol::Real=1e-9`: Tolerance for detecting new points and for Newton iterations.

# Returns
- `Fiber`: A `Fiber` object containing
    - the polynomial `P`,
    - the list of points in the fiber,
    - the basepoint in the codomain, and
    - the critical value bound used for lifting.

# Description
The function computes the fiber over a basepoint `y0` chosen outside the 
critical value disk. A loop around the "point at infinity" in the codomain
is lifted to the domain using monodromy, collecting all preimages of `y0`.
Duplicate points are avoided by a relative tolerance check.

This method guarantees that the returned fiber contains exactly `degree(P)` 
points. Errors are raised if the fiber cannot be fully constructed 
or if any path lifting fails.

# Example
```julia
P = PolynomialCover([1.0, 0.0, -1.0])  # P(z) = z^2 - 1
F = fiber(P; nl=16, s=4, tol=1e-12)
println("Fiber points: ", F)
```
"""
function fiber(
	P::PolynomialCover{T};
	nl::Int = 8,
	s::Real = 3,
	tol::Real = 1e-9,
) where {T<:Real}
	Rcrit = min_critical_radius(coefficients(P))
	x0 = basepoint(P; s=s)

	eS = bound_values_radius_origin(P, Rcrit)
	w0 = value(P, x0)

	loop = basepoint_initial_loop(P; nl=nl, s=s)
	d = degree(P)

	states = PointState{T}[PointState(x0, w0)]

	new_point(z) =
		all(abs(z - st.x) > tol*(1 + abs(st.x)) for st in states)

	i = 1
	while i ≤ length(states) && length(states) < d
		st = states[i]

		st_new = PointState(st.x, st.w)

		st_new, ok = lift_path_fast!(
			P, st_new, loop;
			S=eS,
			newton_tol=tol,
		)

		ok || error("Lifting failed during fiber computation")

		new_point(st_new.x) && push!(states, st_new)

		i += 1
	end

	length(states) == d || error("Fiber not complete")

	xs = [st.x for st in states]

	return Fiber(P, xs, w0, eS)
end

using Base.Threads

function transport_threads!(F::Fiber, path::Path; kwargs...)

    states = [PointState(x, F.w) for x in F.xs]

    Threads.@threads for i in eachindex(states)
        lift_path_fast!(F.P, states[i], path; S=F.S, kwargs...)
    end

    @inbounds for i in eachindex(states)
        F.xs[i] = states[i].x
    end

    F.w = path.vertices[end].z

    return F
end


function transport_serial!(F::Fiber, path::Path; kwargs...)
    # Build states at current base value
    states = [PointState(x, F.w) for x in F.xs]

    # Lift each sheet independently
    for st in states
        lift_path_fast!(F.P, st, path; S=F.S, kwargs...)
    end

    # Update fiber
    F.xs .= (st.x for st in states)
    F.w = path.vertices[end].z

    return F
end

"""
    transport!(F::Fiber, path::Path; parallel=true, kwargs...) -> Fiber

Transport a fiber along a given path in the codomain, updating the fiber points
using path lifting of the polynomial cover.

# Arguments
- `F::Fiber`: The fiber to be transported. Contains the polynomial `P`, the current
  points of the fiber `xs`, the codomain value `w`, and the critical value bound `S`.
- `path::Path`: The path in the codomain along which to lift the fiber.
- `parallel::Bool=true`: Whether to lift each sheet in parallel using Julia threads
  (requires `Threads.nthreads() > 1`).
- `kwargs...`: Additional keyword arguments are forwarded to `lift_path_fast!`,
  e.g., `newton_tol` or `maxiter`.

# Returns
- `Fiber`: The updated fiber after transporting along the path. The `xs` field
  contains the lifted points, and `w` is updated to the final point of the path.

# Description
Each sheet of the fiber (preimage of the basepoint) is lifted independently along
the path using `lift_path_fast!`. If `parallel=true` and multiple threads are
available, the lifting is done concurrently for faster execution. Otherwise, a
serial lifting is performed.

The function ensures that after completion, all sheets of the fiber are updated
and consistent with the final point of the path.

# Example
```julia
using PolynomialCovers

P = PolynomialCover([1.0, 0.0, -1.0])   # P(z) = z^2 - 1
F = fiber(P)
γ = safe_path_arc(basepoint(P), -basepoint(P))

# Transport fiber along the path (automatically selects parallel if available)
transport!(F, γ)

println("Fiber after transport: ", F)

"""

function transport!(F::Fiber, path::Path; parallel=true, kwargs...)
    if parallel && Threads.nthreads() > 1
        transport_threads!(F, path; kwargs...)
    else
        transport_serial!(F, path; kwargs...)
    end
end


"""
    transport_to!(F::Fiber, y::Number; kwargs...) -> Fiber

Transport a fiber to a specified codomain value `y` along a safe path
that avoids the critical values of the polynomial.

# Arguments
- `F::Fiber`: The fiber to transport. Contains the polynomial `P`, current
  points `xs`, current codomain value `w`, and critical value bound `S`.
- `y::Number`: Target value in the codomain to transport the fiber to.
- `kwargs...`: Additional keyword arguments forwarded to `transport!`, e.g.,
  `parallel=true`, `newton_tol`, etc.

# Returns
- `Fiber`: The updated fiber after transport. The `xs` field contains the lifted points,
  and `w` is updated to `y`.

# Description
The function first checks that `y` is sufficiently far from the critical value
disk: |y| must exceed `3S/2` (where `S` is the radius enclosing all critical values).
It then constructs a safe path from the current fiber value `F.w` to `y` using
`safe_path_arc` and calls `transport!` to lift each sheet along this path.

# Errors
Throws an `ArgumentError` if `|y| ≤ 3S/2`, since transport too close to the
critical values may fail or produce incorrect fibers.

# Example
```julia
using PolynomialCovers

P = PolynomialCover([1.0, 0.0, -1.0])
F = fiber(P)
y_target = 3 * critical_values_radius(P)

# Transport fiber to y_target safely
transport_to!(F, y_target)

println("Fiber at y = ", y_target, ": ", F)
"""


function transport_to!(F::Fiber, y::Number; kwargs...)

    if abs(y) ≤ F.S*(3/2)
        throw(ArgumentError(
            "Target y too close to the critical value disk. " *
            "Require |y| > 3S/2 for safe transport."
        ))
    end

    path = safe_path_arc(F.w, y, F.S)
    return transport!(F, path; kwargs...)
end



"""
    cyclic_lift_solver(P::PolynomialCover{T}, y::Number; s::Real=3, tol::Real=1e-9, kwargs...) where {T<:Real}

Compute the cyclically ordered fiber of the polynomial `P` over a target `y` outside the critical value disk.

# Arguments
- `P::PolynomialCover{T}`: the polynomial cover
- `y::Number`: target value in ℂ
- `s::Real=3`: safety factor for the critical value disk
- `tol::Real=1e-9`: tolerance for numerical solvers
- `kwargs...`: additional keyword arguments for `transport_to!`

# Returns
- `Vector{Complex{T}}`: roots of `P(x) = y` in cyclic order around the base fiber

# Notes
- Requires `|y| > 3/2 * S`, where `S` is the radius of the critical value disk (possibly inflated by `s`).
- If `y` is too close to the critical values, an error is raised.
"""

function cyclic_lift_solver(
    P::PolynomialCover{T},
    y::Number;
    s::Real = 3,
    tol::Real = 1e-9,
    kwargs...
) where {T<:Real}

    # Convert to complex
    yc = complex(y)

    # Estimate critical value disk radius
    Rcrit = min_critical_radius(coefficients(P))               # minimal critical radius
    S = real(bound_values_radius_origin(P, Rcrit))             # outer bound on |critical values|
    
    # Safety margin to avoid tangency issues
    margin = 1 - 10eps(real(S))

    # Check if y is safely outside the critical value disk
    if abs(yc) <= 1.5 * S * margin
        error("""
cyclic_lift_solver requires |y| > 3S/2.
Computed S = $S
|y| = $(abs(yc))

Use the general fiber solver instead.
""")
    end

    # Compute the base fiber
    F = fiber(P; s=s, tol=tol)

    # Transport fiber to target yc
    transport_to!(F, yc; kwargs...)

    # Return the roots in cyclic order
    return F.xs
end

#=

function cyclic_lift_solver(
    P::PolynomialCover{T},
    y::Number;
    s::Real = 3,
    tol::Real = 1e-9,
    kwargs...
) where {T<:Real}

    # Ensure complex target value
    yc = complex(y)

    # Critical radii bounds
    Rcrit= min_critical_radius(coefficients(P))
    S = real(bound_values_radius_origin(P, Rcrit))
    # Robust safety check
    margin = 1 - 10eps(real(S))
    #exprintln("Computed S = $S, |y| = $(abs(yc)), margin = $margin")
    if abs(yc) <= 1.5*S * margin
        error("""
cyclic_lift_solver requires |y| > 3S/2.
Computed S = $S
|y| = $(abs(yc))

Use the general fiber solver instead.
""")
    end

    # Base fiber and lifting
    F = fiber(P; s=s, tol=tol)
    transport_to!(F, yc; kwargs...)

    return F.xs   # raíces en orden cíclico
end

=#

const PLOT_DIR = joinpath(@__DIR__, "MyPlots", "PlotExamples")

