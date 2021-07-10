#---------------------------------------------------------------------#
#This function computes the Grid and Elements in 1D.
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function compute_norm(qn,qe,Npoin,DFloat)
    norm2=norm(qn-qe,2)/norm(qe,2)
    return (norm2)
end
