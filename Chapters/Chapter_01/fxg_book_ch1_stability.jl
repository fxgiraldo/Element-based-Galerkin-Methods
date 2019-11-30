using Plots
pyplot()

#---------------------------------------------------------------------#
#This code plots the Stability of for the 1D Wave Equation.
#Written by F.X. Giraldo on 03/2019
#           Department of Applied Maths
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#Variables:
#---------------------------------------------------------------------#
#method=input(" 1=FE & Upwind \n 2=FE & Centered \n 3=FE & Downwind \n 4=BE & Upwind \n 5=BE & Centered \n 6=BE & Downwind \n 7=TR & Upwind \n 8=TR & Centered \n 9=TR & Downwind \n Enter Method: ")
method=1

n=100
dtheta=2*pi/n
theta=0:dtheta:2*pi
if (method <= 3)
    cmax=1
elseif [method >= 4]
    cmax=10
end
m=10*cmax
dc=cmax/m
c=0:dc:cmax
em1=exp.(-im*theta)
ep1=exp.(+im*theta)

a=ones(size(em1))
#Plot Stability of Specific Equations
xmatrix=Array{Float64}(undef,n,m+2)
ymatrix=Array{Float64}(undef,n,m+2)
if method == 1 #FE + Upwind
    for j=1:m+1
        z=a - c[j]*(a - em1) #1st order upwind [explicit 1st order in time]
        xmatrix[:,j]=real(z)
        ymatrix[:,j]=imag(z)
    end
elseif method == 2 #FE + Centered
    for j=2:m+1
        z=a - 0.5*c[j]*(ep1 - em1); #2nd order centered [explicit 1st order in time]
        xmatrix[:,j]=real(z)
        ymatrix[:,j]=imag(z)
    end
elseif method == 3 #FE + Downwind
    for j=2:m+1
        z=a - c[j]*(ep1 - a); #1st order downwind [explicit 1st order in time]
        xmatrix[:,j]=real(z)
        ymatrix[:,j]=imag(z)
    end
elseif method == 4 #BE + Upwind
    for j=2:m+1
        z=a/(a + c[j]*(a - em1)) #1st order uwpind, (implicit 1st order in time)
        xmatrix[:,j]=real(z)
        ymatrix[:,j]=imag(z)
    end
elseif method == 5 #BE + Centered
    for j=2:m+1
        z=a/(a + 0.5*c[j]*(ep1 - em1)) #2nd order centered, (implicit 1st order in time)
        xmatrix[:,j]=real(z)
        ymatrix[:,j]=imag(z)
    end
elseif method == 6 #BE + Downwind
    for j=2:m+1
        z=a/(a + c[j]*(ep1 - a)) #1st order downwind, (implicit 1st order in time)
        xmatrix[:,j]=real(z)
        ymatrix[:,j]=imag(z)
    end
elseif method == 7 #TR + Upwind
    for j=2:m+1
        z=(a - 0.5*c[j]*(a - em1))./(a + 0.5*c[j]*(a - em1)) #1st order upwind [Trapezoidal 2nd order in time]
        xmatrix[:,j]=real(z)
        ymatrix[:,j]=imag(z)
    end
elseif method == 8 #TR + Centered
    for j=2:m+1
        z=(a - 0.25*c[j]*(ep1 - em1))./(a + 0.25*c[j]*(ep1 - em1)) #2nd order centered [TR 2nd order in time]
        xmatrix[:,j]=real(z)
        ymatrix[:,j]=imag(z)
    end
elseif method == 9 #TR + Downwind
    for j=2:m+1
        z=(a - 0.5*c[j]*(ep1 - a))./(a + 0.5*c[j]*(ep1 - a)) #1st order downwind [TR 2nd order in time]
        xmatrix[:,j]=real(z)
        ymatrix[:,j]=imag(z)
    end
end

#Plot Unit Circle
z=ep1
xmatrix[:,m+2]=real(z)
ymatrix[:,m+2]=imag(z)
#plot_handle=plot(xmatrix,ymatrix,leg=false,w=2,aspect_ratio=1)
plot(xmatrix,ymatrix,leg=false,w=1,line=(:dash,1),aspect_ratio=1)
plot!(real(z),imag(z),leg=false,w=1,color=[:black],aspect_ratio=1)
xaxis!("Re(z)")
yaxis!("Im(z)")
title!("Von Neumann Stability Analysis")
#display(plot(plot_handle))
