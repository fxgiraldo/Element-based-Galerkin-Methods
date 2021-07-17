#---------------------------------------------------------------------#
#This function computes the Norm of two arrays
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function compute_norms(qn,qe,DFloat)

    norms=zeros(DFloat,2,1)
    norms[1]=norm(qn[:,1]-qe[:,1],2)/norm(qe[:,1],2)
    norms[2]=norm(qn[:,2]-qe[:,2],2)/norm(qe[:,2],2)
    return (norms)

end
