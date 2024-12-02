#---------------------------------------------------------------------#
#This function computes the global Mass and Differentiation Matrices.
#Written by F.X. Giraldo on April 19 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function create_rhs(qp,M,De,Dwe,intma,periodicity,cgdg_method,flux_method,DFloat)

    include("create_rhs_strong.jl")
    include("create_rhs_weak.jl")

    #Build RHS vector
    if cgdg_method=="strong"
        R=create_rhs_strong(qp,M,De,intma,periodicity,flux_method,DFloat)
    elseif cgdg_method=="weak"
        R=create_rhs_weak(qp,M,Dwe,intma,periodicity,DFloat)
    end

    #Return arrays
    return(R)
end
