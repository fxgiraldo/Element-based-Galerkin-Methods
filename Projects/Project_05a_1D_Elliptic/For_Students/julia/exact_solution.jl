#---------------------------------------------------------------------#
#This function computes the Initial & Exact Solutions.
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function exact_solution(coord,Npoin,case,DFloat)

    #Set some constants
    c=2*π

    #Initialize
    qe=zeros(DFloat,Npoin)
    qeₓ=zeros(DFloat,Npoin)
    fe=zeros(DFloat,Npoin)

    #Generate Grid Points
    for I=1:Npoin
        x=coord[I]
        if (case == 1) #Gaussian
            qe[I]=sin(c*x)
            qeₓ[I]=-c*cos(c*x)
            fe[I]=-c^2*sin(c*x)
        end
    end #I
    return (qe,qeₓ,fe)
end
