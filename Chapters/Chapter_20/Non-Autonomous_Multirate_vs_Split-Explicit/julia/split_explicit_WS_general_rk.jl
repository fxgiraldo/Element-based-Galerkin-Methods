#---------------------------------------------------------------------#
#This routine advances the solution in time using Split-Explicit using
#a simple General Order RK method.
#Written by F.X. Giraldo / P.R. Mugg on 9/21/19
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#
function split_explicit_WS_general_rk(q0,c,I,M,Δt,time)

    #build RK Coefficients
    α𝑖 = zeros(I+1,1)
    for i = 1:I
        α𝑖[i+1] = 1.0/(I-i+1)
    end
    c𝑖 = α𝑖

    #Outer RK3 Loop
    Qi = q0
    Δτ = Δt/M
    for i = 2:I+1
        Qk = q0
        t = time + Δt*c𝑖[i-1]
        (rhs,f𝑖,g𝑖) = rhs_function(Qi,c,t)

        #Inner Euler Loop
        for k = 1:c𝑖[i]*M
            t += Δτ*(k-1)
            (rhs,f𝑘,g𝑘) = rhs_function(Qk,c,t)
            Qk += Δτ*(f𝑘 + g𝑖)
        end
        Qi = Qk
    end

    return qp = Qi
end
