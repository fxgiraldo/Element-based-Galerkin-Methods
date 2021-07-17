#=
---------------------------------------------------------------------
This function computes the Metric Terms for General 2D Quad Grids.

The strategy follows Algorithm 12.3 in F.X. Giraldo's Introduction to Element-based 
Galerkin Methods using Tensor-Product Bases: Analysis, Algorithms, and Applications.

Written by F.X. Giraldo on 5/2019
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

include("map_deriv.jl")

function compute_metrics(coord,intma,ψ,dψ,Ne,Np,Nq,DFloat)

    #Initialize Global Arrays
    ξ_x=zeros(DFloat,Nq,Nq,Ne)
    ξ_y=zeros(DFloat,Nq,Nq,Ne)
    η_x=zeros(DFloat,Nq,Nq,Ne)
    η_y=zeros(DFloat,Nq,Nq,Ne)
    jac=zeros(DFloat,Nq,Nq,Ne)

    #Initialize Local Arrays
    x=zeros(DFloat,Np,Np)
    y=zeros(DFloat,Np,Np)

    #loop thru the elements
    for e=1:Ne

        #Store Element Variables
        for j=1:Np
            for i=1:Np
                I=intma[i,j,e]
                x[i,j]=coord[1,I]
                y[i,j]=coord[2,I]
            end #i
        end #j

        #Construct Mapping Derivatives: dx/dξ; dx/dη; dy/dξ; dy/dη
        (x_ξ,x_η)=map_deriv(ψ,dψ,x,Np,Nq,DFloat)
        (y_ξ,y_η)=map_deriv(ψ,dψ,y,Np,Nq,DFloat)

        #Construct Inverse Mapping: dξ/dx; dξ/dy; dη/dx; dη/dy
        for j=1:Nq
            for i=1:Nq
                xjac=x_ξ[i,j]*y_η[i,j] - x_η[i,j]*y_ξ[i,j]
                ξ_x[i,j,e]=+1.0/xjac*y_η[i,j]
                ξ_y[i,j,e]=-1.0/xjac*x_η[i,j]
                η_x[i,j,e]=-1.0/xjac*y_ξ[i,j]
                η_y[i,j,e]=+1.0/xjac*x_ξ[i,j]
                jac[i,j,e]=abs(xjac)
            end #i
        end #j
    end #e

    return (ξ_x,ξ_y,η_x,η_y,jac)
end
