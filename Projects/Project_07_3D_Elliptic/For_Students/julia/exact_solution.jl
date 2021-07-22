#---------------------------------------------------------------------#
#This function computes the Initial & Analytic Solutions.
#Written by F.X. Giraldo on 7/2021
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function exact_solution(coord,Npoin,c,DFloat)

    #Initialize
    qe=zeros(DFloat,Npoin,1)
    fe=zeros(DFloat,Npoin,1)
    cc=c*π

    #Generate Grid Points
    for I=1:Npoin
        x=coord[1,I]
        y=coord[2,I]
        z=coord[3,I]
        qe[I]=sin(cc*x)*sin(cc*y)*sin(cc*z)
        fe[I]=-3*cc^2*sin(cc*x)*sin(cc*y)*sin(cc*z)
    end #I

    return (qe,fe)

end


