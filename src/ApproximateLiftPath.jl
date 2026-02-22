########################################
# ApproximateLiftPath.jl
########################################




#Path(trace::Trace{T}) where {T<:Real} = Path(trace.vertices)

#=
function lift_segment_linear!(
    P::PolynomialCover,
    state::LiftState{T},
    w1::Complex{T},
    w2::Complex{T},
) where {T<:Real}

    Δw = w2 - w1
    dP = derivative(P, state.x)
    @assert abs(dP) > 0 "Derivative vanished outside critical disk"

    x_new = state.x + Δw / dP

    state.x = x_new
    state.w = value(P, x_new)

    return nothing
end


function aux_lift_path!(
    P::PolynomialCover,
    state::LiftState{T},
    ypath::Path{T};
    S::Real,
) where {T<:Real}

    trace = state.xpath
    i = length(trace.vertices)

    # No podemos avanzar más de lo que hay en ypath
    if i ≥ length(ypath.vertices)
        return state, ypath, true
    end

    x_i   = state.x
    y_i   = ypath.vertices[i].z
    y_ip1 = ypath.vertices[i+1].z

    dP = derivative(P, x_i)
    if abs(dP) == 0
        println("DEBUG aux_lift_path!: x_i=", x_i, " y_i=", y_i, " y_ip1=", y_ip1, " dP=", dP)
        println("DEBUG: P.coeffs=", P.coeffs, " degree=", degree(P), " typeof(x_i)=", typeof(x_i), " eltype(coeffs)=", eltype(P.coeffs))
        error("Derivative vanished outside critical disk")
    end

    x_new = x_i + (y_ip1 - state.w) / dP
    w_new = value(P, x_new)

    # Caso final: estamos intentando cerrar el camino
    if i + 1 == length(ypath.vertices)
        if !segment_intersects_disk(w_new, y_ip1, 3*S/2)
            push!(trace.vertices, Vertex(x_new))
            state.x = x_new
            state.w = w_new
            return state, ypath, true
        else
            y_mid = (y_i + y_ip1) / 2
            insert!(ypath.vertices, i+1, Vertex(y_mid))
            return state, ypath, false
        end
    end

    # Caso general
    y_ip2 = ypath.vertices[i+2].z

    safe1 = !segment_intersects_disk(w_new, y_ip1, 3*S/2)
    safe2 = !segment_intersects_disk(w_new, y_ip2, 3*S/2)

    if safe1 && safe2
        push!(trace.vertices, Vertex(x_new))
        state.x = x_new
        state.w = w_new
        return state, ypath, false
    else
        y_mid = (y_i + y_ip1) / 2
        insert!(ypath.vertices, i+1, Vertex(y_mid))
        return state, ypath, false
    end
end


function refine_endpoint_newton!(
    P::PolynomialCover,
    state::LiftState{T},
    w_target::Complex{T};
    tol::Real = 1e-9,
    max_iter::Int = 100,
) where {T<:Real}

    x = state.x

    for _ in 1:max_iter
        Px = value(P, x)
        err = Px - w_target

        if abs(err) ≤ tol
            state.x = x
            state.w = Px
            return x
        end

        dPx = derivative(P, x)
        abs(dPx) == 0 && error("Derivative vanished during Newton refinement")

        x -= err / dPx
    end

    error("Newton refinement did not converge in $max_iter iterations")
end

=#

function aux_lift_predictor!(
    P::PolynomialCover,
    state::PointState{T},
    ypath::Path{T};
    S::Real,
    max_subdiv::Int = 10,
) where {T<:Real}

    x = state.x
    w = state.w

    verts = ypath.vertices
    n = length(verts)

    i = 1
    subdiv = 0

    while i < n

        y_i   = verts[i].z
        y_ip1 = verts[i+1].z

        # predictor step
        dP = derivative(P, x)
        abs(dP) < eps(T) && error("Near critical point")

        x_new = x + (y_ip1 - w)/dP
        w_new = value(P, x_new)

        # homotopy safety check
        if segment_intersects_disk(w_new, y_ip1, 3S/2)

            subdiv += 1
            subdiv > max_subdiv && return state, false

            # halve step in codomain
            y_mid = (y_i + y_ip1)/2

            # do predictor only for half-step
            x = x + (y_mid - w)/dP
            w = value(P, x)

            continue
        end

        # accept step
        x = x_new
        w = w_new
        i += 1
    end

    state.x = x
    state.w = w

    return state, true
end

function refine_endpoint_newton!(
    P::PolynomialCover,
    state::PointState{T},
    w_target::Complex{T};
    tol::Real = 1e-9,
    max_iter::Int = 100,
) where {T<:Real}

    x = state.x

    for _ in 1:max_iter
        Px = value(P, x)
        err = Px - w_target

        if abs(err) ≤ tol*(1 + abs(w_target))
            state.x = x
            state.w = w_target
            return x
        end

        dPx = derivative(P, x)
        abs(dPx) < sqrt(eps(T)) &&
            error("Derivative nearly vanished during Newton refinement")

        x -= err/dPx
    end

    error("Newton refinement did not converge in $max_iter iterations")
end

#=
function lift_path_with_trace!(
    P::PolynomialCover,
    state::LiftState{T},
    ypath::Path{T};
    S::Real,
    newton_tol::Real = 1e-12,
) where {T<:Real}

    finished = false
    while !finished
        _, ypath, finished = aux_lift_path!(P, state, ypath; S=S)
    end

    xfinal = refine_endpoint_newton!(
        P,
        state,
        ypath.vertices[end].z;
        tol = newton_tol,
    )

    push!(state.xpath.vertices, Vertex(xfinal))

    # Convertimos la traza en Path *solo aquí*
    xpath = Path(state.xpath.vertices)

    return (
        xpath = xpath,
        ypath = ypath,
        xfinal = xfinal,
    )
end
=#
function lift_path_fast!(
    P::PolynomialCover,
    state::PointState{T},
    ypath::Path{T};
    S::Real,
    newton_tol::Real = 1e-9,
) where {T<:Real}

    # --- predictor lifting (homotopy-safe) ---
    state_pred, ok = aux_lift_predictor!(P, state, ypath; S=S)

    ok || return state_pred, false

    # --- endpoint Newton refinement ---
    w_target = ypath.vertices[end].z
    refine_endpoint_newton!(
        P,
        state_pred,
        w_target;
        tol=newton_tol,
    )

    return state_pred, true
end


