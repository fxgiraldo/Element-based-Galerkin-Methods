#---------------------------------------------------------------------#
#This function computes the time-step size for a given max Courant number
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
"""
    compute_courant(q0,qb,Courant_max,dx_min,DFloat)

returns: dt and courant
"""
function compute_courant(q0,Courant_max,dx_min,DFloat)

    #Recompute DT
    Npoin=size(q0,1)
    umax=DFloat(0.0)
    for i=1:Npoin
        u=q0[i]
        umax=max(umax, abs(u))
    end
    dt=Courant_max*dx_min/umax
    courant=umax*dt/dx_min
    #println("Courant = ",courant," dt = ",dt)
    return (dt,courant)
end
