using Test
using PolynomialCovers


#using Polynomials
#using PolynomialRoots 

function same_fiber(A, B; atol=1e-8, rtol=1e-8)
    length(A) == length(B) || return false
    used = falses(length(B))

    for a in A
        found = false
        for (j, b) in enumerate(B)
            if !used[j] && isapprox(a, b; atol=atol, rtol=rtol)
                used[j] = true
                found = true
                break
            end
        end
        found || return false
    end
    return true
end# o quita esto si usas roots de Polynomials.jl


P = PolynomialCover([1.0, 1.0, 1.0, 1.0])
Q= PolynomialCover([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0])

# 1. radio crítico mínimo (geométrico)
    coeff_P=coefficients(P)
    println("coeff_P: ", coeff_P)
    Rcrit= PolynomialCovers.min_critical_radius(coeff_P)
    Rcrit_Q=PolynomialCovers.min_critical_radius(coefficients(Q))

    # 2. punto base en el dominio (seguro)
    x0 = PolynomialCovers.basepoint(P)
    x0_Q =PolynomialCovers. basepoint(Q)

    # 3. cota en valores
    S = PolynomialCovers.bound_values_radius_origin(P, Rcrit)
    S_Q = PolynomialCovers.bound_values_radius_origin(Q, Rcrit_Q)
    Swork_Q = 3S_Q/2
    Swork = 3S/2

    # 4. valor base
    w0 = value(P, x0)
    w0_Q = value(Q, x0_Q)

    # 5. lazo seguro en el codominio
    loop1 = PolynomialCovers.basepoint_initial_loop(P)
    loop1_Q = PolynomialCovers. basepoint_initial_loop(Q)

    d = degree(P)
    d_Q = degree(Q)


    pointstate1 = PolynomialCovers.PointState(x0, w0)
    pointstate1_Q = PolynomialCovers.PointState(x0_Q, w0_Q) 
  
    fvc=PolynomialCovers.fiber_via_companion(P, w0)
    

    lfap=PolynomialCovers.aux_lift_predictor!(P, pointstate1, loop1; S = S)
    lfap_Q=PolynomialCovers.aux_lift_predictor!(Q, pointstate1_Q, loop1_Q; S = S_Q)
 


lf1=PolynomialCovers.lift_path_fast!(
    P,
    pointstate1,
    loop1; S = S, newton_tol = 1e-9)
lf1_Q=PolynomialCovers.lift_path_fast!(
    Q,
    pointstate1_Q,
    loop1_Q; S = S_Q, newton_tol = 1e-9)
    
    typeof(P)
    typeof(Q)


    

    fiber_cover =PolynomialCovers.fiber(P)

    fiber_cover.P
    fiber_cover.xs
 
       
    fiber_cover.w
    fiber_cover.S

    fiber_cover_Q =  fiber_cover =PolynomialCovers.fiber(Q)
    fiber_cover_Q.P
    fiber_cover_Q.xs
    fiber_cover_Q.w
    fiber_cover_Q.S
    
    
  
    # polinomio CORRECTO
    fr = PolynomialCovers.fiber_via_companion(P, w0)
    fr_Q = PolynomialCovers.fiber_via_companion(Q, w0_Q)
    fc = fiber_cover.xs

    fiber_cover =PolynomialCovers.cyclic_lift_solver(P, 0.0+60.0im )


@testset "PolynomialCover: fiber vs roots" begin
    P = PolynomialCover([1.0, 2.0, 3.0, 4.0])

    fiber_cover = PolynomialCovers.fiber(P)
    fc = fiber_cover.xs

    #@test !isempty(fc)

    # valor base REAL del fiber
    w = value(P, fc[1])

    # polinomio CORRECTO
  
    fr = PolynomialCovers.fiber_via_companion(P, w)

    @test same_fiber(fc, fr; atol=1e-8, rtol=1e-8)
end

@testset "PolynomialCover: fiber vs roots" begin
    Q= PolynomialCover([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0])

    fiber_cover = PolynomialCovers.fiber(Q)
        fc = fiber_cover.xs

    #@test !isempty(fc)

    # valor base REAL del fiber
    w = value(Q, fc[1])

    # polinomio CORRECTO
  
    fr = PolynomialCovers.fiber_via_companion(Q, w)

    @test same_fiber(fc, fr; atol=1e-8, rtol=1e-8)
end