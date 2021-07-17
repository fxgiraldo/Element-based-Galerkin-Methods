#---------------------------------------------------------------------#
#This function computes the Norm of two vectors
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
using LinearAlgebra

function compute_norms(qn,qe,npoin,DFloat)

    norms=zeros(DFloat,3)
    norms[1]=norm(qn-qe,1)/norm(qe,1)
    norms[2]=norm(qn-qe,2)/norm(qe,2)
    norms[3]=norm(qn-qe,Inf)/norm(qe,Inf)
    return (norms)

end
