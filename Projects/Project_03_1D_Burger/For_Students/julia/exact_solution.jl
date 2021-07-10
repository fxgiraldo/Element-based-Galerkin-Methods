#---------------------------------------------------------------------#
#This function computes the Initial & Exact Solutions.
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function exact_solution(coord,Npoin,time,case,DFloat)

    #Set some constants
    xmin=0
    xmax=+2
    xc=0.5*(xmax+xmin)
    xl=xmax-xmin
    rc=0.125
    sigma=1/16; #best for right convergence rates
    # sigma=1/8
    # sigma=1/4; #fatter Gaussian to compare against 2D code
    
    #Initialize
    qe=zeros(DFloat,Npoin)

    #timec=time - floor(time)
    timec=time
    #Generate Grid Points
    for I=1:Npoin
        x=coord[I]
        if (case == 1) #Sin
            qe[I]=sin(pi*x) + 0.01;
        end
    end #I
    return (qe)
end
