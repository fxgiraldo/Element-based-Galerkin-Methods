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
    plot_handle=plot(coord,q,xlabel="x",ylabel="q(x,t)",legend=true,lw=3,label=[space_method],title=text)
    display(plot_handle)

end
