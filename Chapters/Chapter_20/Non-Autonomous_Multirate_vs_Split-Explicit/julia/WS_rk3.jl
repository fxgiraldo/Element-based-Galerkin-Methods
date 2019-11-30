#---------------------------------------------------------------------#
#This routine advances the solution in time using Split-Explicit using
#a simple General Order RK method.
#Written by F.X. Giraldo / P.R. Mugg on 9/20/19
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#
function WS_rk3(q0,c,Δt,time)

    #1st RK Stage
    Q = q0
    (rhs,f,g) = rhs_function(Q,c,time)
    Q = q0 + Δt/3*(f + g)

    #2nd RK Stage
    (rhs,f,g) = rhs_function(Q,c,time + Δt/3)
    Q = q0 + Δt/2*(f + g)

    #3rd RK Stage
    (rhs,f,g) = rhs_function(Q,c,time + Δt/2)
    Q = q0 + Δt*(f + g)

    #Update
    return qp = Q

end #end function
