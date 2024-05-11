#---------------------------------------------------------------------#
#This function computes the Initial & Analytic Solutions.
#Written by F.X. Giraldo on 10/2003
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function exact_solution(coord,Npoin,time,icase,DFloat)

    #Initialize
    qe=zeros(DFloat,Npoin,1)
    ue=zeros(DFloat,2,Npoin)
    
    #Constants
    xmin=minimum(coord[1,:]) 
    xmax=maximum(coord[1,:])
    ymin=minimum(coord[2,:]) 
    ymax=maximum(coord[2,:])
    xl=xmax-xmin
    yl=ymax-ymin
    xm=0.5*(xmax+xmin)
    ym=0.5*(ymax+ymin)
    xc=xmin + 0.25*xl
    yc=ymin + 0.5*yl
    σ=32.0
   
    #Initial Condition
    for I=1:Npoin
        x=coord[1,I]
        y=coord[2,I]

        if (icase == 1)  #Gaussian in CCW direction
            ue[1,I]=+(y-ym);
            ue[2,I]=-(x-xm);
            xx=x - xc*cos(time) - yc*sin(time)
            yy=y + xc*sin(time) - yc*cos(time)
            qe[I]=exp( - σ*(xx^2 + yy^2 ) ) 
        end
    end #I

    return (qe,ue)

end


