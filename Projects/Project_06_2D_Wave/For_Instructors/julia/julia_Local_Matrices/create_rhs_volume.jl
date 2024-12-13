#=
---------------------------------------------------------------------
This function computes the Volume integral contribution to the RHS vector
for a unified CGDG method for the 2D Wave Equation.

Written by F.X. Giraldo on July 12, 2021
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

function create_rhs_volume(q,u,De_x,De_y,intma,periodicity,DFloat)
 
    #Get lengths of arrays
    (Np,Np,Ne)=size(intma)
    Npoin=length(q)

    #Store local arrays
    rhs=zeros(DFloat,Npoin)

    #build differentiation matrix contribution
    for e=1:Ne
        for i=1:Np, j=1:Np
            I=intma[i,j,e]
            for k=1:Np, l=1:Np
                K=intma[k,l,e]
                rhs[periodicity[I]]+=De_x[i,j,k,l,e]*q[K]*u[1,K] + De_y[i,j,k,l,e]*q[K]*u[2,K]
            end #k,l
        end #i,j
    end #e
    
    #Return arrays
    return(rhs)
end
