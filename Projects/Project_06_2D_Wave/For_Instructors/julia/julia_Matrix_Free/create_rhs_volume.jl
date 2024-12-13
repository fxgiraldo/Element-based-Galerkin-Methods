#=
---------------------------------------------------------------------
This function computes the Volume integral contribution to the RHS vector
for the Weak Form unified CGDG method for the 2D Wave Equation.

Written by F.X. Giraldo on December 13, 2024
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

function create_rhs_volume(q,u,ψ,dψ,ξ_x,ξ_y,η_x,η_y,jac,ωq,intma,periodicity,DFloat)
 
    #Get lengths of arrays
    (Np,Np,Ne)=size(intma)
    Npoin=length(q)

    #Store local arrays
    rhs=zeros(DFloat,Npoin)

    #build differentiation matrix contribution
    for e=1:Ne
        for i=1:Np, j=1:Np
            I=intma[i,j,e]
            weight=ωq[i]*ωq[j]*jac[i,j,e]
            e_x=ξ_x[i,j,e]
            e_y=ξ_y[i,j,e]
            n_x=η_x[i,j,e]
            n_y=η_y[i,j,e]

            #Interpolate at Integration Point
            I=intma[i,j,e]
            u_k=u[1,I]
            v_k=u[2,I]
            q_k=q[I]
            f_k=q_k*u_k
            g_k=q_k*v_k

            for k=1:Np
                #derivative of flux along ξ
                I=intma[k,j,e]
                dhqde=dψ[k,i]*ψ[j,j]*( e_x*f_k + e_y*g_k )
                rhs[periodicity[I]] += weight*dhqde
                
                #derivative of flux along η
                I=intma[i,k,e]
                dhqdn=ψ[i,i]*dψ[k,j]*( n_x*f_k + n_y*g_k )
                rhs[periodicity[I]] += weight*dhqdn
            end #k
        end #i,j
    end #e
    
    #Return arrays
    return(rhs)
end
