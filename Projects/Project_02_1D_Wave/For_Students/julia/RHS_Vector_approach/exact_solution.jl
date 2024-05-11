#---------------------------------------------------------------------#
#This function computes the Initial & Exact Solutions.
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function exact_solution(coord,Npoin,time,case,DFloat)

    #Set some constants
    xmin=-1
    xmax=+1
    xc=0.5*(xmax+xmin)
    xl=xmax-xmin
    rc=0.125
    sigma=1/16; #best for right convergence rates
    # sigma=1/8
    # sigma=1/4; #fatter Gaussian to compare against 2D code
    u=DFloat(2.0)

    #Initialize
    qe=zeros(DFloat,Npoin)

    #timec=time - floor(time)
    timec=time
    #Generate Grid Points
    for I=1:Npoin
        x=coord[I]
        xbar=xc + u*timec
        if (xbar > xmax)
            xbar=xmin + (xbar-xmax)
        end
        r=x-xbar
        if (case == 1) #Gaussian
            qe[I]=exp( -(x-xbar)^2/(2*sigma)^2 )
        elseif (case == 2) #Square Wave
      	    if ( abs(r) <= rc )
                qe[I]=1
            end
        end
    end #I
    return (qe,u)
end
