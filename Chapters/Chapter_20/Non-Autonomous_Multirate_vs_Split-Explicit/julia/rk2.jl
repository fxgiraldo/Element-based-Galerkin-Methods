#-----------------------------------------------------------------#
# This routine advances the solution in time using Midpoint RK2
# Written by F.X. Giraldo / P.R. Mugg on 9/19/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#------------------------------------------------------------------#

function rk2(q0,c,Δt,time)

    qi = q0

    #Evaluate at Beginning of Interval
    (rhs,f,g) = rhs_function(q0,c,time)

    #update to Midpoint of Interval
    qi = q0 + (Δt/2)*rhs

    #update to End of Interval
    (rhs,f,g) = rhs_function(qi,c,time+(Δt/2))
    qi = q0 + Δt*rhs

    #update
    return qp = qi
end
