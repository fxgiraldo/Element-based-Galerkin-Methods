#---------------------------------------------------------------------#
#This function performs the time-integration
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
include("ti_SSP.jl")
include("ti_LSRK.jl")

function time_integration!(q0,u,coord,Dhat,Fhat,intma,periodicity,time,ntime,dt,stages,DFloat,ti_type,space_method,plot_movie)

    if (ti_type == "SSP")
        q0 = ti_SSP!(q0,u,coord,Dhat,Fhat,intma,periodicity,time,ntime,dt,stages,space_method,plot_movie,DFloat)
    elseif (ti_type == "LSRK")
        q0 = ti_LSRK!(q0,u,coord,Dhat,Fhat,intma,periodicity,time,ntime,dt,space_method,plot_movie,DFloat)
    end
    return (q0)
end
