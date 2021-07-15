#=
---------------------------------------------------------------------
This function computes the normal vectors and Jacobian of each FACE.

Written by F.X. Giraldo on July 14, 2021
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

include("map_deriv.jl")

function compute_normals(face,intma,coord,Nface,Np,Nq,ψ,dψ)

    #global arrays
    normals=zeros(DFloat,2,Nq,Nface)
    jac_face=zeros(DFloat,Nq,Nface)

    #local arrays
    x=zeros(DFloat,Np,Np)
    y=zeros(DFloat,Np,Np)

    #loop thru the sides
    for f=1:Nface

        #Store Left & Right Elements
        ilocl=face[1,f]
        ilocr=face[2,f]
        iel=face[3,f]
        ier=face[4,f]

        #Store Element Variables
        for j=1:Np, i=1:Np
            I=intma[i,j,iel]
            x[i,j]=coord[1,I]
            y[i,j]=coord[2,I]
        end #i, j
 
        #Construct Mapping Derivatives: dx/dksi; dx/deta; dy/dksi
        #dy/deta
        (x_ξ,x_η)=map_deriv(ψ,dψ,x,Np,Nq,DFloat)
        (y_ξ,y_η)=map_deriv(ψ,dψ,y,Np,Nq,DFloat)

        #Compute Normals 
        for l=1:Nq
            if (ilocl == 1) 
                #Side 1: eta=-1
                i=l
                j=1
                normals[1,l,f]=+y_ξ[i,j]
                normals[2,l,f]=-x_ξ[i,j]
                jac_face[l,f]=sqrt( x_ξ[i,j]^2 + y_ξ[i,j]^2 )
            elseif (ilocl == 2) 
                #Side 2: ksi=+1
                i=Nq
                j=l
                normals[1,l,f]=+y_η[i,j]
                normals[2,l,f]=-x_η[i,j]
                jac_face[l,f]=sqrt( x_η[i,j]^2 + y_η[i,j]^2 )
            elseif (ilocl == 3)
                #Side 3: eta=+1
                i=Nq+1-l
                j=Nq
                normals[1,l,f]=-y_ξ[i,j]
                normals[2,l,f]=+x_ξ[i,j]
                jac_face[l,f]=sqrt( x_ξ[i,j]^2 + y_ξ[i,j]^2 )
            elseif (ilocl == 4)
                #Side 4: ksi=-1
                i=1
                j=Nq+1-l
                normals[1,l,f]=-y_η[i,j]
                normals[2,l,f]=+x_η[i,j]
                jac_face[l,f]=sqrt( x_η[i,j]^2 + y_η[i,j]^2 )
            end
        end #l 

        #Normalize Normals
        for l=1:Nq
            rnx=sqrt( normals[1,l,f]^2 + normals[2,l,f]^2 )
            normals[1,l,f]=normals[1,l,f]/rnx
            normals[2,l,f]=normals[2,l,f]/rnx
        end #l   

    end #s
    
    return(normals,jac_face)
end
