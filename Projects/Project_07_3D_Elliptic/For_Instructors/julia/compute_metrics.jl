#=
---------------------------------------------------------------------
This function computes the Metric Terms for General 3D Hexahedral Grids.

The strategy follows Algorithm 12.4 and 12.5 in F.X. Giraldo's Introduction to Element-based 
Galerkin Methods using Tensor-Product Bases: Analysis, Algorithms, and Applications.
Here, we use the cross-product metric terms. As an exercise the student may be asked to build the 
Curl-Invariant metric terms.

Written by F.X. Giraldo on 7/2021
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

include("map_deriv.jl")

function compute_metrics(coord,intma,ψ,dψ,Ne,Np,Nq,DFloat)

    #Initialize Global Arrays
    ξ_x=zeros(DFloat,Nq,Nq,Nq,Ne)
    ξ_y=zeros(DFloat,Nq,Nq,Nq,Ne)
    ξ_z=zeros(DFloat,Nq,Nq,Nq,Ne)
    η_x=zeros(DFloat,Nq,Nq,Nq,Ne)
    η_y=zeros(DFloat,Nq,Nq,Nq,Ne)
    η_z=zeros(DFloat,Nq,Nq,Nq,Ne)
    ζ_x=zeros(DFloat,Nq,Nq,Nq,Ne)
    ζ_y=zeros(DFloat,Nq,Nq,Nq,Ne)
    ζ_z=zeros(DFloat,Nq,Nq,Nq,Ne)
    jac=zeros(DFloat,Nq,Nq,Nq,Ne)

    #Initialize Local Arrays
    x=zeros(DFloat,Np,Np,Np)
    y=zeros(DFloat,Np,Np,Np)
    z=zeros(DFloat,Np,Np,Np)

    #loop thru the elements
    for e=1:Ne

        #Store Element Variables
        for k=1:Np, j=1:Np, i=1:Np
            I=intma[i,j,k,e]
            x[i,j,k]=coord[1,I]
            y[i,j,k]=coord[2,I]
            z[i,j,k]=coord[3,I]
        end #j

        #Construct Mapping Derivatives: dx/dξ; dx/dη; dy/dξ; dy/dη
        (x_ξ,x_η,x_ζ)=map_deriv(ψ,dψ,x,Np,Nq,DFloat)
        (y_ξ,y_η,y_ζ)=map_deriv(ψ,dψ,y,Np,Nq,DFloat)
        (z_ξ,z_η,z_ζ)=map_deriv(ψ,dψ,z,Np,Nq,DFloat)

        #Construct Inverse Mapping: dξ/dx; dξ/dy; dη/dx; dη/dy
        for k=1:Nq, j=1:Nq, i=1:Nq
            xj = (x_ξ[i,j,k]*y_η[i,j,k]*z_ζ[i,j,k] - x_ξ[i,j,k]*y_ζ[i,j,k]*z_η[i,j,k]) -
            (y_ξ[i,j,k]*x_η[i,j,k]*z_ζ[i,j,k] - y_ξ[i,j,k]*x_ζ[i,j,k]*z_η[i,j,k]) +
            (z_ξ[i,j,k]*x_η[i,j,k]*y_ζ[i,j,k] - z_ξ[i,j,k]*x_ζ[i,j,k]*y_η[i,j,k])  

            ξ_x[i,j,k,e]= (y_η[i,j,k]*z_ζ[i,j,k]-y_ζ[i,j,k]*z_η[i,j,k])/xj
            ξ_y[i,j,k,e]=-(x_η[i,j,k]*z_ζ[i,j,k]-x_ζ[i,j,k]*z_η[i,j,k])/xj
            ξ_z[i,j,k,e]= (x_η[i,j,k]*y_ζ[i,j,k]-x_ζ[i,j,k]*y_η[i,j,k])/xj
            η_x[i,j,k,e]=-(y_ξ[i,j,k]*z_ζ[i,j,k]-y_ζ[i,j,k]*z_ξ[i,j,k])/xj
            η_y[i,j,k,e]= (x_ξ[i,j,k]*z_ζ[i,j,k]-x_ζ[i,j,k]*z_ξ[i,j,k])/xj
            η_z[i,j,k,e]= -(x_ξ[i,j,k]*y_ζ[i,j,k]-x_ζ[i,j,k]*y_ξ[i,j,k])/xj
            ζ_x[i,j,k,e]= (y_ξ[i,j,k]*z_η[i,j,k] -y_η[i,j,k]*z_ξ[i,j,k])/xj
            ζ_y[i,j,k,e]=-(x_ξ[i,j,k]*z_η[i,j,k] -x_η[i,j,k]*z_ξ[i,j,k])/xj
            ζ_z[i,j,k,e]= (x_ξ[i,j,k]*y_η[i,j,k] -x_η[i,j,k]*y_ξ[i,j,k])/xj   
            jac[i,j,k,e]=abs(xj)      
        end #i, j, k
    end #e

    return (ξ_x,ξ_y,ξ_z,η_x,η_y,η_z,ζ_x,ζ_y,ζ_z,jac)
end
