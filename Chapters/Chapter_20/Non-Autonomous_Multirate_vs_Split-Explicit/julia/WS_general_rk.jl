#---------------------------------------------------------------------#
#This routine advances the solution in time using
#a simple General Order RK method.
#Written by F.X. Giraldo / P.R. Mugg on 9/21/19
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#
function WS_general_rk(q0,c,I,Δt,time)

    #build RK Coefficients
    αi = zeros(I+1,1)
    for i = 1:I
        αi[i+1] = 1.0/(I-i+1)
    end

    qi = q0
    (rhs,f,g) = rhs_function(qi,c,time)
    for i = 1:I
        #Update
        qi = q0 + Δt*αi[i+1]*rhs

        #Construct RHS
        (rhs,f,g) = rhs_function(qi,c,(time + Δt*αi[i+1]))
    end #for i

    #final update
    return qp = qi

end #end function
