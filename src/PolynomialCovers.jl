module PolynomialCovers

include("PreLifting.jl")



# Estructura principal para la fibra
mutable struct Fiber{T<:Real}
	P::PolynomialCover{T}
	xs::Vector{Complex{T}}
	w::Complex{T}
	S::Real
end

# Construcción de la fibra usando monodromía en el infinito
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


function transport!(F::Fiber, path::Path; kwargs...)
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
    transport_to!(F, y)

Safely transport the fiber to a new value `y`,
avoiding the critical value disk.
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
    cyclic_lift_solver(P, y; s=3, tol=1e-9, kwargs...)

Compute the roots of P(x)=y using fast monodromy lifting from a base fiber.

This solver assumes that `|y| > 3S/2`, where `S` bounds the critical
value disk of the polynomial cover. Under this condition the lifting
is fast and numerically stable.

Returns the roots in cyclic order.
"""
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



export PolynomialCover, coefficients, value, degree,
	Fiber, fiber, transport!, transport_to!, cyclic_lift_solver

end # module
