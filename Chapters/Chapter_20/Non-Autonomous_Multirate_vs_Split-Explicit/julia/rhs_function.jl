#-----------------------------------------------------------------#
# This routine computes the RHS forcing of a simple two-rate ODE
# Written by F.X. Giraldo / P.R. Mugg on 9/19/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#------------------------------------------------------------------#
function rhs_function(q0,c,time)

    # RHS components
    f = c*q0;        #fast component
    g = sin(time);   #slow component

    # Construct RHS
    rhs = f + g;

    return (rhs,f,g)
end #end of rhs_function
