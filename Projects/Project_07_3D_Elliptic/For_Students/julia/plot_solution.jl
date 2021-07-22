#---------------------------------------------------------------------#
#This function plots the solution
#Written by F.X. Giraldo on April 19 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function plot_solution(q,coord,space_method,time)

    #Plot Solution
    text=string(space_method,", time= ","$time")
    p1=plot(coord,q[1,:],xlabel="x",ylabel="h(x,t)",legend=true,lw=3,label=[space_method],title=text,ylims=(-0.5,0.5))
    p2=plot(coord,q[2,:],xlabel="x",ylabel="U(x,t)",legend=true,lw=3,label=[space_method],title=text,ylims=(-0.5,0.5))
    display(plot(p1,p2))
    #display(p1)

end
