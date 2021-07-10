#---------------------------------------------------------------------#
#This function plots the solution
#Written by F.X. Giraldo on April 19 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function plot_solution(q0,qb,coord,space_method,time,case)

    #Exact Solution
    Npoin=size(q0,1)
    (qe,qb,gravity) = exact_solution(coord,Npoin,time,case,DFloat)
    
    #Plot Solution
    p1 = plot(coord,qe[:,1],linecolor = :blue,label = "Analytic Solution",legend = :bottomright,legendfont=font(5))
    plot!(p1,coord,q0[:,1],linestyle = :dashdot,linecolor = :red,label = "$space_method Solution")
    title!(p1,"Height (hs)")
    xlabel!(p1,"x")
    ylabel!(p1,"hs(x,t)")
    plot!(p1,xlims = (first(coord),last(coord)+.1),ylims = (-1,0.5))

    p2 = plot(coord,qe[:,2],linecolor = :blue,label = "Analytic Solution",legend = :bottomright,legendfont=font(5))
    plot!(p2,coord,q0[:,2],linestyle = :dashdot,linecolor = :red,label = "$space_method Solution")
    title!(p2,"Momentum (U)")
    xlabel!(p2,"x")
    ylabel!(p2,"U(x,t)")
    plot!(p2,xlims = (first(coord),last(coord)+.1),ylims = (-1,0.5))

    figure_1 = plot(p1,p2,layout = (1,2))
    display(figure_1)
end
