#=
---------------------------------------------------------------------
This function computes the RHS vector for a unified CGDG method for the 
2D Wave Equation.

Written by F.X. Giraldo on July 12, 2021
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

include("create_rhs_volume.jl")
include("create_rhs_flux.jl")

function create_rhs(q,u,M,De_x,De_y,intma,periodicity,face,mapL,mapR,normals,jac_face,ωq,space_method,DFloat)
 
    #--------------------------------------------------------------------------------------------------------
    #Students Add their Contributions Here: They should construct one function for the Volume Integral 
    #contribution and a separate one for the Flux Integral contribution. This way they can optimize each one 
    #separately
    #--------------------------------------------------------------------------------------------------------

    #Construct Volume integral contribution
    rhs = create_rhs_volume(q,u,De_x,De_y,intma,periodicity,DFloat)

    #Construct Flux integral contribution (will also work for CG but cancels completely due to periodicity)
    if (space_method == "DG")
        rhs = create_rhs_flux!(rhs,q,u,face,normals,jac_face,ωq,mapL,mapR,intma,periodicity,DFloat)
    end

    #Multiply by inverse Mass matrix
    rhs=M\rhs

    #Apply Periodicity
    Npoin=length(q)
    for i=1:Npoin
        rhs[i]=rhs[periodicity[i]]
    end
    
    #Return arrays
    return(rhs)
end
