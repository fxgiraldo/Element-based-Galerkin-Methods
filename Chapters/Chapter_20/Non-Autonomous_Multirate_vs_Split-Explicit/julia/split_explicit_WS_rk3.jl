#---------------------------------------------------------------------#
#This routine advances the solution in time using Split-Explicit using
#a simple General Order RK method.
#Note that for M=3 in 1st stage, M=2 in 2nd, and M=1 in 3rd gives
#3rd order convergence as expected.
#
#As is, we get 0th order convergence but this works in NUMA.
#
#Written by F.X. Giraldo / P.R. Mugg on 9/21/19
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#
function split_explicit_WS_rk3(q0,c,M,Δt,time)

    #1st RK Stage
    Qi = q0
    t = time
    (rhs,f𝑖,g𝑖) = rhs_function(Qi,c,t)

    #M=3
    Δτ = Δt/M
    Qk = q0
    for k = 1:M/3
        t += Δτ*(k-1)
        (rhs,f𝑘,g𝑘) = rhs_function(Qk,c,t)
        #Qk += Δτ*(f𝑘 + g𝑖)
        Qk = q0 + Δτ*(f𝑘 + g𝑖)
    end
    Qi = Qk

    #2nd RK Stage
    #M=2
    t = time + Δt/3
    (rhs,f𝑖,g𝑖) = rhs_function(Qi,c,t)
    dτ = Δt/M
    #Qk = q0
    for k = 1:M/2
        t += dτ*(k-1)
        (rhs,f𝑘,g𝑘) = rhs_function(Qk,c,t)
        #Qk += Δτ*(f𝑘 + g𝑖)
        Qk = q0 + dτ*(f𝑘 + g𝑖)
    end
    Qi = Qk

    #3rd RK Stage
    #M=1
    t = time + Δt/2
    (rhs,f𝑖,g𝑖) = rhs_function(Qi,c,t)
    Δτ = Δt/M
    #Qk = q0
    for k = 1:M
        t += Δτ*(k-1)
        (rhs,f𝑘,g𝑘) = rhs_function(Qk,c,t)
        #Qk += Δτ*(f𝑘 + g𝑖)
        Qk = q0 + Δτ*(f𝑘 + g𝑖)
    end

    return qp = Qk
end #end function
