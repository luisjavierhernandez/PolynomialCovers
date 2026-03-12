module PolynomialCovers


# =========================
# Basic structures
# =========================

include("Geometry.jl")

include("PolynomialUtils.jl")

include("PolynomialCriticalDisks.jl")

include("PreLifting.jl")

include("ApproximateLiftPath.jl")


include("FiberTransport.jl")



#========================                    
export PolynomialCover, coefficients, value, degree, Path,
 PointState, Disk, basepoint, basepoint_initial_loop,
	Fiber, fiber,PointState, transport!, transport_to!, cyclic_lift_solver,
	Vertex, safe_path_arc, critical_radius, 
    y_sequence, aux_lift_predictor!, segment_intersects_disk,
    min_critical_radius, bound_values_radius_origin, 
    critical_values_radius_origin,  critical_values_disk_origin,    
     critical_values_exclusion_radius, _min_radius, avoid_points_path,
     startpoint,endpoint,PLOT_DIR, plot_path

=======================#                    


export PolynomialCover,
       coefficients,
       degree,
       value,

       Vertex,
       Path,
       Disk,
       distance,

       Fiber,
       critical_values_radius_origin,
       critical_values_disk_origin,
       fiber,
       PointState,

       basepoint,
       basepoint_initial_loop,

       safe_path_arc,
       avoid_points_path,

       transport!,
       transport_to!,
       cyclic_lift_solver
     

end # module
