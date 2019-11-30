#---------------------------------------------------------------------#
#This function computes the Grid and Elements in 1D.
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function compute_norm(qn,qe,Npoin,DFloat)

    norm=zeros(DFloat,1)
    l2_top=zeros(DFloat,1)
    l2_bot=zeros(DFloat,1)
    for i=1:Npoin
        l2_top[1]=l2_top[1] + (qn[i]-qe[i])^2
        l2_bot[1]=l2_bot[1] + qe[i]^2
    end
    #    norm[1]=sqrt( l2_top[1] /l2_bot[1] )
    norm[1]=sqrt( l2_top[1] )

    return (norm)
end
